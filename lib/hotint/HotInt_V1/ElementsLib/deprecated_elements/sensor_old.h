//#**************************************************************
//#
//# filename:             sensor_old.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
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

#ifndef SENSOR__H
#define SENSOR__H

#include "sensors.h"

////////////////////////////////////////////////////////
// This is an old and obsolete sensor implementation. //
////////////////////////////////////////////////////////

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++       MBSSensor          +++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef enum {TSempty = 0, TSElement = 1, TSNode = 2, TSPos = 4, TSVel = 8, TSAngle = 16,
TSX = 32, TSY = 64, TSZ = 128, TSDist = 256, TSDeflection = 512, TSDOF = 1024, TSplanar = 2048, TSAuxElem = 4096,
TSAccel = 8192, /*TSMultSensor = 16384, */TSLocalAxis = 32768, TSForce = 65536, TSStress = 131072, TSStrain = 262144, 
TSInelStrain = 524288, TSOutputSensor = 1048576, TSCurvature = 2097152, TSOtherLocalAxisElement = 4194304, 
TSKardanAngle = 8388608, TSDisplacement = 16777216, TSEigenValue = 33554432, TSConstraintDrift = 67108864, TSFFT = 134217728, /*TSActorMoment = (use TSForce!!!)*/
TSSpecialSensorValue = 268435456, TSAxialExtension = 536870912, TSXData = 1073741824, TSLoad = 2147483648, TSFile = 4294967296} TMBSSensor;

//possible sensors:
//TSDOF ... degree of freedom sensor
//TSElement + {TSPos|TSVel|TSAccel} + {TSX|TSY|TSZ} //measure x/y/z component of position/velocity/acceleration; for this sensor, a local position for measuring must be added
//TSElement + {TSAngle} + {TSX|TSY|TSZ} {+TSVel} {+ TSLocalAxis} //measure angle w.r.t. X/Y/Z (local) axis; for this sensor, a local position for measuring must be added
																																//TSVel: measure angular velocity
//TSElement + TSplanar + {TSPos|TSVel|TSAccel} + {TSX|TSY}                     //measure x/y/z component of position/velocity; for this sensor, a local position for measuring must be added
//TSElement + TSplanar +  TSAngle {+TSVel} //measure angle axis; for this sensor, a local position for measuring must be added; TSVel: measure angular velocity
//TSElement + {TSKardanAngle} + {TSX|TSY|TSZ}  //measure cardan angle w.r.t. TSX=A(gamma,z),TSY=A(gamma,z).A(beta,y) or TSZ=A(gamma,z).A(beta,y).A(alpha,x) rotation axe; for this sensor, the cardan angles are computed by the euler parameters of a BODY3D
//TSFFT ...  additionally flag: create file "SENSORFILENAME-fft.txt" with columns "frequency | amplitude | phase", (set 'ifdef 1' at all comments '// RL: in test phase' to try fft computations) 
//TSFile... read t-y, values from file, use function SetTSFile for definition of the file name, columns etc.
typedef enum {TSMempty = 0, TSMAverage = 1, TSMSum = 2, TSMMult = 4,
TSMDiv = 8, TSMNorm = 16, TSMDifference = 32, TSMMax =64, TSMMin = 128} TMBSMultSensor; //norm: x1^2+x2^2+x3^2 ...; Difference is done with factor=-1 !!!

typedef enum {TSMises = 0, TSCompXX = 1, TSCompYY = 2, TSCompZZ = 3, 
TSCompYZ = 4, TSCompXZ = 5, TSCompXY = 6} TMBSTensorcomponent;

typedef enum 
{
	TSCnothing = 0,		//do nothing
	TSCmin = 1,				//output minimum
	TSCmax = 2,				//output maximum
	TSCminmax = 3,		//output min and max
	TSCamplitude = 4, //output max-min
	TSCmaxabs = 5,    //maximum of fabs(min) and fabs(max)
	TSCmeanL2 = 6,    //mean L2-norm
} TSensorComputation;


//$ RL 2011-7-14: FFT code, 
// these flags are used in combination to the flags TSensorComputation and define further options in case of fft based computation value 
typedef enum 
{
	TSCFFTnothing = 0,          //do nothing
	TSCFFTamplitude = 1,        //analyze amplitude spectrum for evaluation of computation value
	TSCFFTphase = 2,            //analyze phase spectrum  for evaluation computation value (ignored in case of TSCmeanL2)
	TSCFFTlocalExtremeVals = 4, //local max./min. (searched from left to right) instead of global max./min. (ignored in case of  TSCmeanL2); use SetNumberOfLocalExtremeValues for activation of this flag and for definition of number of extreme value(s)
	TSCFFTfrequency = 8         //output frequency based cost function value (e.g. frequency at (local) max./min); (ignored, if TSCmeanL2 is active)
} TSensorComputationFFT;


// as an intermediate stage before the final transition to the new sensor,
// the old MBSSensor is now subclassed from the new base Sensor class;
class MBSSensor : public Sensor
{
public:
	MBSSensor(): elementlist(), nodelist(), poslist(), poslist2D(), fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(NULL)
	{
		Init();
	}
	MBSSensor(MBS* mbsi, TMBSSensor typeI): elementlist(), nodelist(), poslist(), poslist2D(), fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
	}
	//element number and local node number
	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element, int node): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		elementlist.Add(element);
		nodelist.Add(node);

		SetSensorName(GetTypeName());
	}
	//node number: global node number; for pos/vel use TSPos/TSVel; for direction use TSX, TSY,TSZ
	//only tested for FiniteElement3D, only for displacement
	// if type = TSEigenValue --> "node" is index of global eigenvalue vector
	MBSSensor(MBS* mbsi, TMBSSensor typeI, int node): elementlist(), nodelist(), poslist(), poslist2D(), fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		nodelist.Add(node);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element, const Vector3D& pos, TMBSTensorcomponent tensor_compI = TSMises): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		tensor_comp = tensor_compI;
		elementlist.Add(element);
		poslist.Add(pos);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element, const Vector3D& pos1, const Vector3D& pos2): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		//e.g. for TSAngle:
		Init();
		type = typeI;
		elementlist.Add(element);
		poslist.Add(pos1); 
		poslist.Add(pos2);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element, const Vector3D& pos1, const Vector3D& pos2, const Vector3D& pos3, const Vector3D& pos4=Vector3D(0.)): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		//e.g. for TSAngle:
		Init();
		type = typeI;
		elementlist.Add(element);
		poslist.Add(pos1); 
		poslist.Add(pos2);
		poslist.Add(pos3);
		poslist.Add(pos4);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element1, int element2, const Vector3D& pos1, const Vector3D& pos2, const Vector3D& pos3=Vector3D(0.), 
		const Vector3D& pos4=Vector3D(0.), const Vector3D& pos5=Vector3D(0.), const Vector3D& pos6=Vector3D(0.)): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		//e.g. for TSAngle:
		Init();
		type = typeI;
		elementlist.Add(element1);
		elementlist.Add(element2);
		poslist.Add(pos1); 
		poslist.Add(pos2);
		poslist.Add(pos3);
		poslist.Add(pos4);
		poslist.Add(pos5);
		poslist.Add(pos6);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, const TArray<int>& elements, const TArray<int>& nodes): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		elementlist = elements;
		nodelist = nodes;

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, const TArray<int>& elements, const TArray<Vector3D>& positions): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		elementlist = elements;
		poslist = positions;

		SetSensorName(GetTypeName());
	}

	// D.R. 12.01.2011 e.g. for TSDeflection
	MBSSensor(MBS* mbsi, TMBSSensor typeI, const TArray<int>& elements, const TArray<int>& nodes, const TArray<Vector3D>& positions): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		elementlist = elements;
		nodelist = nodes;
		poslist = positions;

		SetSensorName(GetTypeName());
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//2D:
	MBSSensor(MBS* mbsi, TMBSSensor typeI, int element, const Vector2D& pos, TMBSTensorcomponent tensor_compI = TSMises): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		tensor_comp = tensor_compI;

		elementlist.Add(element);
		poslist2D.Add(pos);

		SetSensorName(GetTypeName());
	}

	MBSSensor(MBS* mbsi, TMBSSensor typeI, const TArray<int>& elements, const TArray<Vector2D>& positions): elementlist(), nodelist(), poslist(), poslist2D(),
		fft_frequencies(), fft_amplitudes(), fft_phase(), Sensor(mbsi)
	{
		Init();
		type = typeI;
		elementlist = elements;
		poslist2D = positions;

		SetSensorName(GetTypeName());
	}

	//------------------------------------
	// old: don't use this function anymore!!!
	void SetSensorComputation(int sensorComputationFlag, double startTimeSensorComputationI = -1e30, double endTimeSensorComputationI = 1e30) // in case of no parameter, the max/min-values are computed the whole simulation
	{
		GetMBS()->UO(UO_LVL_warn) << mystr("WARNING: type of SensorComputation parameter 1 changed from int to TSensorComputation, please change in model file\n");
		
		SetSensorComputation( (TSensorComputation) sensorComputationFlag, startTimeSensorComputationI, endTimeSensorComputationI);
	}
	//------------------------------------

	
	void SetSensorComputation(TSensorComputation typeI, double startTimeSensorComputationI = -1e30, double endTimeSensorComputationI = 1e30) // in case of no parameter, the max/min-values are computed the whole simulation
	{
		startTimeSensorComputation = startTimeSensorComputationI;
		endTimeSensorComputation = endTimeSensorComputationI;
		
		minVal =  1.0e30;
		maxVal = -1.0e30;
		refVal = 0;
		n_sensor_computations = 0;

		SetFlag_SensorComputation(typeI);		
	}

	//$ RL 2011-7-14: FFT code, 
	// set parameters for costfunction evaluations in frequency domain
	void SetSensorComputationFFT(TSensorComputation sctypeI, TSensorComputationFFT fft_flagI, double fft_constant_sample_timeI, double startFreqencySensorComputationI = -1e30, double endFrequencySensorComputationI = 1e30, double startTimeSensorComputationI = -1e30, double endTimeSensorComputationI = 1e30, int NLocalExtremeValues = 0) // in case of no parameter, the max/min-values are computed the whole simulation
	{						
		if((NLocalExtremeValues && (sctypeI & TSCmeanL2 || sctypeI & TSCmaxabs || sctypeI & TSCamplitude)) || (!NLocalExtremeValues && fft_flagI & TSCFFTlocalExtremeVals))
		{
			GetMBS()->UO(UO_LVL_err) << "ERROR: local extreme values not defined.";
			return;
		}
		
		if(fft_flagI & TSCFFTlocalExtremeVals && sctypeI & TSCmeanL2)
		{
			GetMBS()->UO(UO_LVL_err) << "ERROR: local extreme values not defined.";
			return;
		}

		if(fft_flagI & TSCFFTphase && sctypeI & TSCmeanL2)
		{
			GetMBS()->UO(UO_LVL_err) << "ERROR: L2-norm for phase spectrum not defined.";
			return;
		}
		if(!(type & TSFFT))
		{
			type = (TMBSSensor)(type + TSFFT);
			GetMBS()->UO(UO_LVL_warn) << "WARNING: Sensor attribute TSFFT added.";
		}

		if(fft_constant_sample_timeI<=0.)
		{
			fft_constant_sample_time = 1e-3;
			GetMBS()->UO(UO_LVL_warn) << "WARNING: FFT-sample time changed to 0.001 s.";
		}
		else
		{
			fft_constant_sample_time = fft_constant_sample_timeI;
		}
	
		SetSensorComputation(sctypeI, startTimeSensorComputationI, endTimeSensorComputationI);
		SetFlag_SensorComputationFFT(fft_flagI);
		ModifySignalStorageMode(SSM_InternalArray);
		startFrequencySensorComputation = startFreqencySensorComputationI; 
		endFrequencySensorComputation = endFrequencySensorComputationI;

		if(NLocalExtremeValues)
		{			
			SetNumberOfLocalExtremeValues(NLocalExtremeValues);
		}
	}

	MBSSensor(const MBSSensor& s) : Sensor(NULL)
	{
		CopyFrom(s);
	};

	/*MBSSensor& operator=(const MBSSensor& s) 
	{
		if (this == &s) {return *this;}
		CopyFrom(s);
		return *this;
	}*/

	virtual MBSSensor* GetCopy()
	{
		MBSSensor* ec = new MBSSensor();
		ec->CopyFrom(*this);
		return ec;
	}

	virtual void CopyFrom(const MBSSensor& s)
	{
		Init();

		Sensor::CopyFrom(s);

		type = s.type;
		elementlist = s.elementlist;
		nodelist = s.nodelist;
		poslist = s.poslist;
		poslist2D = s.poslist2D;

		offset = s.offset;
		factor = s.factor;
		
		tensor_comp = s.tensor_comp;

		// calculate min/max value during simulation
		minVal = s.minVal;
		maxVal = s.maxVal;
		refVal = s.refVal;
		n_sensor_computations = s.n_sensor_computations;
		startTimeSensorComputation = s.startTimeSensorComputation;
		endTimeSensorComputation = s.endTimeSensorComputation;
//$ RL 2011-7-14: FFT code
		startFrequencySensorComputation = s.startFrequencySensorComputation;
		endFrequencySensorComputation = s.endFrequencySensorComputation;



		flag_sensor_computation = s.flag_sensor_computation;
		computationValueWeight = s.computationValueWeight;
		// file sensor values

		if(s.tsfileValues)
		{		
			tsfileValues = s.tsfileValues->GetCopy();
		}
		tsfileName = s.tsfileName;                // name of data file for file sensor TSFile

		// reference values
		hasReferenceValues = s.hasReferenceValues;
		if(hasReferenceValues && s.referenceValues)
		{
			referenceValues = s.referenceValues->GetCopy();
		}
		referenceSensorNumber = s.referenceSensorNumber;

		//times = s.times;
		//values = s.values;
//$ RL 2011-7-14: FFT code		
		fft_fileout = s.fft_fileout;
		fft_constant_sample_time = s.fft_constant_sample_time;
		flag_sensor_computation_fft = s.flag_sensor_computation_fft;
		fft_NLocalExtremeValues = s.fft_NLocalExtremeValues;
		fft_comp_diff_of_ampl = s.fft_comp_diff_of_ampl;
		fft_frequencies = s.fft_frequencies;
		fft_amplitudes = s.fft_amplitudes;
		fft_phase = s.fft_phase;
	}

	virtual ~MBSSensor() {
		if(referenceValues){delete referenceValues; referenceValues = 0;};
		delete tsfileValues; tsfileValues = 0;
		/*if(fileout){fileout->close();delete fileout;fileout=0;}*/  //$ YV 2012-06: will be closed by mbs
		if(fft_fileout){fft_fileout->close();delete fft_fileout;fft_fileout=0;}
	}

// initialize sensor before constructor enters data
	virtual void Init()
	{
		type = TSempty;
		tensor_comp = TSMises;
		offset = 0;
		factor = 1;

		actvalue = 0;
		
		hasReferenceValues = 0;
		referenceValues = 0;
		// calculate min/max value during simulation
		double inf = 1e300;
		SetMinVal(inf); //there is always at least one value
		SetMaxVal(-inf); //there is always at least one value
		SetRefVal(0.);
		n_sensor_computations = 0;
		SetStartTimeSensorComputation(-inf);
		SetEndTimeSensorComputation(inf);
		fft_frequencies.SetLen(0);
		fft_amplitudes.SetLen(0);
		fft_phase.SetLen(0);

//$ RL 2011-7-14: FFT code
		SetStartFrequencySensorComputation(-inf);
		SetEndFrequencySensorComputation(inf);




		SetFlag_SensorComputation(TSCnothing);
		tsfileValues = 0;
		tsfileName = mystr("");                // name of data file for file sensor TSFile

		SetReferenceSensor(0); // initially, no reference sensor added
		//TArray<double> weights;
		//weights.Add(1.);
		//SetSensorComputationValueWeights(weights);
		SetSensorComputationValueWeight(1.);


//$ RL 2011-7-14: FFT code
		fft_fileout = 0;	
		fft_constant_sample_time = 0.001; // default value
		SetFlag_SensorComputationFFT(TSCFFTnothing);
		fft_NLocalExtremeValues = 0;
		fft_comp_diff_of_ampl = 0;
	}

	virtual bool IsConsistent(mystr& errorstr); //$ YV 2012-06: the old logic with three states was not used anyway -> //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer	//$ DR 2012-10 return value changed from bool to int

	virtual mystr GetTypeName(); //get name for the corresponding type
	virtual void SetReferenceSolution(const MathFunction& refval) 
	{
		referenceValues = new MathFunction(refval);
		hasReferenceValues = 1;
	}
	virtual int HasReferenceValues() const {return hasReferenceValues;}                        // are reference values subtracted from actvalue
	virtual void GetReferenceValues(const Vector& times, Vector& refvalues) const;             // reference data vector from MathFunction are evaluated at given time-points "times" and stored in the vector "refvalues"

  //---------------------------------------------------------------------------------
	// set reference sensor for computation of difference
	virtual void SetReferenceSensor(const int referenceSensorNumberi) //$ RL 2011-02:[ 
	{		
		referenceSensorNumber = referenceSensorNumberi;
		if(referenceSensorNumber)
		{
			GetMBS()->GetSensor(referenceSensorNumber).ModifySignalStorageMode(SSM_InternalArray); // store time and value data in local storage of MBSSensor (fft computations possible)
		}
	}                                                                 //$ RL 2011-02:]
	virtual int GetReferenceSensorNumber() const{return referenceSensorNumber;}
  //---------------------------------------------------------------------------------

	/*
	virtual int NElements() const {return elementlist.Length();}
	virtual int GetElNum(int i) const {return elementlist(i);}
	virtual int& GetElNum(int i) {return elementlist(i);}
	*/

	virtual int GetNumberOfRelatedElements() { return elementlist.Length(); }
	virtual int& GetRelatedElementNumber(int nElement) { return elementlist(nElement); }

	virtual void SetFactor(double fact) {factor = fact;}
	virtual double GetFactor() const {return factor;}
	virtual void SetOffset(double off) {offset = off;}
	virtual double GetOffset() const {return offset;}

	//$ YV 2012-06: the sensors may produce just one scalar value; the evaluation procedure looks differently
	virtual double GetCurrentValue(double time); 

	virtual int GetNumberOfDrawingPositions();
	virtual Vector3D GetDrawPosition(int i);
	// virtual void Draw();		//$ YV 2012-06: drawing will be performed in the base class

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // sensor computation functions - e.g. computing maximal and minimal values
	virtual void SetMinVal(double valI) {minVal = valI;}
	virtual double GetMinVal() const {return minVal;}
	virtual void SetMaxVal(double valI) {maxVal = valI;} 
	virtual double GetMaxVal() const {return maxVal;}

	//$ RL 2011-7-14: FFT code, 
	virtual Vector GetExtremeValuesFFT(const Vector& vec, TArray<int>& indices, int minmax = 1) const;  // search one or N extreme values; vec...data vector., ind...indices of extreme values, minmax... 1.0 --> search maxima, -1.0 --> search minima


	virtual void SetRefVal(double valI) {refVal = valI;} 
	virtual double GetRefVal() const {return refVal;}
	virtual Vector GetSensorComputationValue() const;

	virtual void SetSensorComputationValueWeight(double weight){computationValueWeight.SetLen(0);computationValueWeight.Add(weight);} //set cost function weight (1D-Sensor)
	virtual void SetSensorComputationValueWeights(TArray<double>& weight){computationValueWeight.CopyFrom(weight);}                   //set cost function weights
	
	virtual void SetStartTimeSensorComputation(double valI) {startTimeSensorComputation = valI;}
	virtual double GetStartTimeSensorComputation() const {return startTimeSensorComputation;}
	virtual void SetEndTimeSensorComputation(double valI){endTimeSensorComputation = valI;	}
	virtual double GetEndTimeSensorComputation() const {return endTimeSensorComputation;}

	//$ RL 2011-7-14: FFT code
	virtual void SetStartFrequencySensorComputation(double valI) {startFrequencySensorComputation = valI;}
	virtual double GetStartFrequencySensorComputation() const {return startFrequencySensorComputation;}
	virtual void SetEndFrequencySensorComputation(double valI){endFrequencySensorComputation = valI;	}
	virtual double GetEndFrequencySensorComputation() const {return endFrequencySensorComputation;}
	virtual double GetConstantSampleTimeFFT() const {return fft_constant_sample_time;}
	virtual void SetConstantSampleTimeFFT(double valI) {type = (TMBSSensor)(type | TSFFT); fft_constant_sample_time=valI;ModifySignalStorageMode(SSM_InternalArray);}


	virtual void SetFlag_SensorComputation(TSensorComputation flag) {flag_sensor_computation = flag;}
	virtual TSensorComputation GetFlag_SensorComputation() const {return flag_sensor_computation;}
//$ RL 2011-7-14: FFT code
	virtual void SetFlag_SensorComputationFFT(TSensorComputationFFT flag)
	{
		if(flag && 
		 (
			((flag & TSCFFTamplitude) && (flag & TSCFFTphase))	||
			(!(flag & TSCFFTamplitude) && !(flag & TSCFFTphase))
		  )
		 )
		{
				GetMBS()->UO(UO_LVL_err).InstantMessageText("Error: Use amplitude XOR phase spectrum for costfunction computation in frequency domain.");
				return;
		}
		flag_sensor_computation_fft = flag;
	}
	virtual TSensorComputationFFT GetFlag_SensorComputationFFT() const {return flag_sensor_computation_fft;}



	virtual void DoSensorComputations();
	virtual void IncrementNSensorComputation() {n_sensor_computations++;}
	//check if time is in range of sensor computation
	virtual int IsSensorComputation(){return GetMBS()->GetTime() >= GetStartTimeSensorComputation() && GetMBS()->GetTime() <=	GetEndTimeSensorComputation();}
	virtual TSensorComputation FlagSensorComputation() const	{	return flag_sensor_computation;	}
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//$ RL 2011-7-14: FFT code
  // load sensor data from file
	//$ YV 2012-06: MBSSensor::LoadSensorData() was used only in FFT for the actual sensor and for the reference one. But both sensors are explicitly switched into SSM_InternalArray mode. Therefore we don’t need the function and can take the values from the internal arrays.
	//virtual int LoadSensorData(TArray<double>& timesi, TArray<double>& valuesi); // load sensor time and sample values, t ... time points, values ... output of sensor
	//old: virtual void LoadSensorDataFFT(TArray<double>& frequencies, TArray<double>& amplitudes, TArray<double>& phase) const; // load sensor fft-data at the end of computation from file, assumption: single sensor fft-outputfile with data frequency, amplitude and phase must exists
	virtual int IsSensorComputationFFT(){return type & TSFFT;} // returns, if the fft-computation flag is set
	virtual void SetOutputFileFFT(ofstream* file) {fft_fileout = file;} // set pointer to the ofstream of the fft-outputfile
	virtual ofstream* GetOutputFileFFT() {return fft_fileout;}          // get pointer to the ofstream of the fft-outputfile
	virtual void SetNumberOfLocalExtremeValues(int val)                 // sets the number of extreme values (maxima or minima), which should be computed
	{
		fft_NLocalExtremeValues = val;
		flag_sensor_computation_fft = (TSensorComputationFFT)(flag_sensor_computation_fft|TSCFFTlocalExtremeVals);
	}
	virtual void SetFlag_Diff_FFTs_MesRef(int activate_difference_of_ffts){fft_comp_diff_of_ampl=activate_difference_of_ffts;}                        // activates computation of difference fft(ymes)-fft(yref) as result of fft-computation
	virtual int GetFlag_Diff_FFTs_MesRef(){return fft_comp_diff_of_ampl;}                        // activates computation of difference fft(ymes)-fft(yref) as result of fft-computation

	virtual int GetNumberOfLocalExtremeValues() const {return fft_NLocalExtremeValues;} // return number of local extreme values
	virtual int SetTSFile(mystr filename, int col1, int col2, int interp = 1, mystr comment = mystr("%")) //set the file name for definition of the date measured with file sensor (flag TSFile), file name inclusive path, col1...time-data, col2 ... y-data
	{
		tsfileName = filename;
		tsfileValues = new MathFunction();
		if(tsfileValues->SetPiecewiseFromFile2(filename, col1, col2, interp, comment))
		{
			delete tsfileValues;
			tsfileValues = 0;
			GetMBS()->UO().InstantMessageText(mystr("Error during TSFile: File " + filename + "not found!"));
			return 1; // error;
		}
		return 0; // ok
	}
	// virtual void StoreSensorTimeValueData(){ times.Add(GetMBS()->GetTime()); values.Add(GetValue()); } // store time point and value of sensor   //$ YV 2012-06: is now in the base class
	virtual Vector& FFT_frequencies(){return fft_frequencies;}// Get reference to fft storage for frequencies (Hz)
	virtual Vector& FFT_amplitudes(){return fft_amplitudes;}	 // Get reference to fft storage for amplitudes
	virtual Vector& FFT_phase(){return fft_phase;}						 // Get reference to fft storage for phase (rad)
	virtual const Vector& FFT_frequencies() const {return fft_frequencies;}// Get reference to fft storage for frequencies (Hz)
	virtual const Vector& FFT_amplitudes() const {return fft_amplitudes;}	 // Get reference to fft storage for amplitudes
	virtual const Vector& FFT_phase() const {return fft_phase;}						 // Get reference to fft storage for phase (rad)
	virtual void PostComputationOperationsFFT(); //do evaluations of fast fourier transformation after computation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual void SetWriteResults(int mode=3) 
	{
		SetSignalStorageMode((SensorSignalStorageMode)mode);
	}
	virtual void SetFlag_InitPlotTool(int flag=1)
	{
		SetOpenSensorWatchPlotToolAtStartUpFlag(flag);
	}
	virtual void SetDrawingDim(Vector3D drawDimension)
	{
		SetDrawDimension(drawDimension);
	}


protected:
	// MBS* mbs;		//$ YV 2012-06: is now in the base class
	TMBSSensor type;
	TMBSTensorcomponent tensor_comp;
	// mystr sensorname;		//$ YV 2012-06: is now in the base class
	mutable TArray<int> elementlist;
	TArray<int> nodelist;
	TArray<Vector3D> poslist;
	TArray<Vector2D> poslist2D;
	// int precision;	//$ YV 2012-06: is now in the base class
	// ofstream* fileout; //only set if writeresults==2 or 3	//$ YV 2012-06: is now in the base class

	// int sensornumber;		//$ YV 2012-06: not used any longer

	double factor;
	double offset;

	// int writeresults;		//$ YV 2012-06: is governed differently by the base class
	// int visible;		//$ YV 2012-06: is now in the base class
	// Vector3D draw_dim;		//$ YV 2012-06: is now in the base class

	//$ YV 2012-06: the sensors may produce just one scalar value
	mutable double actvalue;		// actually measured value of the sensor -- intrinsic for the MBSSensor class, in the base class the value is saved into other variables

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// sensor compuations
	//$ YV 2012-06: the variables below and the corresponding functionality are left "as is", with the new sensors this functionality will be carried by sensor processors
	double minVal, maxVal;	// calculate min/max value during simulation (eg. flag_sensor_computation€{1,2,...,5})
	double refVal; //reference value computed by sensor (e.g. flag_sensor_computation=6: L2-norm)
	double startTimeSensorComputation;  // start time of sensor computation (used for time OR frequency domain computations)
	double endTimeSensorComputation;    // end time of sensor computation (used for time OR frequency domain computations)
//$ RL 2011-7-14: FFT code
	double startFrequencySensorComputation;  // start frequency of sensor computation (used ONLY for frequency domain computations)
	double endFrequencySensorComputation;    // end frequency of sensor computation   (used ONLY for frequency domain computations)	

	int n_sensor_computations; //number of sensor computations, e.g. needed for norm
	TSensorComputation flag_sensor_computation; //0=do nothing, 1=output minimum, 2=output maximum, 3=output min and max, 4=output amplitude, 5=maximum of fabs(min) and fabs(max), 6=mean L2-norm
	TArray<double> computationValueWeight; // multiplicative weight of sensor computation value (standard: 1.0)
	MathFunction* tsfileValues;          // look - up table with values from file containing t-y values, TSFile
	mystr tsfileName;                // name of data file for file sensor TSFile
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// compute difference between actualvalues and referenceValues
	MathFunction* referenceValues;  // look - up table with time reference value
	int hasReferenceValues;        // flag, if difference should be evaluated
  // --------------------------------
	int referenceSensorNumber;     // use sensor for getting reference values //$ RL 2011-02: 
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//$ RL 2011-7-14: FFT code:[ 
	// TArray<double> times;             // storage for time points of sensor		//$ YV 2012-06: is now in the base class
	// TArray<double> values;            // storage values of sensor						//$ YV 2012-06: is now in the base class
	//$ YV 2012-06: fft processing is left "as is", in the future it should be implemented with the new sensors using sensor processors
	Vector fft_frequencies;   // storage for FFT: frequencies (Hz) 
	Vector fft_amplitudes;		// storage for FFT: amplitudes
	Vector fft_phase;         // storage for FFT: phase (rad)
	ofstream* fft_fileout;            // file for fft-results	
	double fft_constant_sample_time;  // constant sample time for fft-computation 
	int fft_NLocalExtremeValues;      // set number of (local) extreme value(s) for cost function computation (sorted from lower to upper frequency)
	TSensorComputationFFT flag_sensor_computation_fft; //0=do nothing, 1=use amplitude for cost function, 2=use phase for cost function, 4
	int fft_comp_diff_of_ampl; //0=compute fft of time signal in sensor file: fft(yfile) (default), 1=if reference signal exists: yfile=ymes-yref, following is computed: fft(yfile+yref)-fft(yref) e.g. to suppress influence of time delays between ymes and yref.
	//$ RL 2011-7-14: FFT code:] 
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//$ AD 2011-08 for PlotTool
	// mystr displayname;		//$ YV 2012-06: does not seem to be actually used
	// int flag_initialize_plottool_at_assemble;	//$ YV 2012-06: is now in the base class
	// int flag_save_final;	//$ YV 2012-06: does not seem to be actually used

	//$ YV 2012-09: the functions below were added for compatibility with the new Sensor class
	virtual bool HasSensorProcessingEvaluationData()
	{
		return flag_sensor_computation || IsSensorComputationFFT();
	}
	virtual bool NeedsFileForPostComputationProcessing()
	{
		return IsSensorComputationFFT();
	}
	TArray<double> signalProcessingEvaluationData;
	virtual TArray<double> & GetSignalProcessingEvaluationData()
	{
		Vector v = GetSensorComputationValue();
		// now we need to convert to a member variable, because we are returning a reference
		signalProcessingEvaluationData.SetLen(v.Length());
		for(int i = 1; i <= v.Length(); i++)
			signalProcessingEvaluationData(i) = v(i);
		return signalProcessingEvaluationData;
	}
	virtual void ApplyPostComputationSensorProcessing(ofstream * outputFile)
	{
		SetOutputFileFFT(outputFile);
		PostComputationOperationsFFT();
		SetOutputFileFFT(NULL);
	}
	void Evaluate(double time)
	{
		Sensor::Evaluate(time);
		DoSensorComputations();
	}
};


//sensor that consists of a list of sensors:
class MBSMultipleSensor: public Sensor
{
public:
	MBSMultipleSensor() : Sensor(NULL)
	{
		offset = 0;
		factor = 1;
	}
	//generate sensor with list of sensors:
	MBSMultipleSensor(MBS* mbsi, TMBSMultSensor multtypeI, TArray<int>& sensorsI): Sensor(mbsi)
	{
		multtype = multtypeI;
		sensors = sensorsI;
		offset = 0;
		factor = 1;

		SetSensorName("MultipleSensor");
	}

	MBSMultipleSensor(const MBSMultipleSensor& s) : Sensor(NULL)
	{
		CopyFrom(s);
	};

	MBSMultipleSensor& operator=(const MBSMultipleSensor& s) 
	{
		if (this == &s) {return *this;}
		CopyFrom(s);
		return *this;
	}
	virtual ~MBSMultipleSensor() {};

	virtual MBSMultipleSensor* GetCopy()
	{
		MBSMultipleSensor* ec = new MBSMultipleSensor();
		ec->CopyFrom(*this);
		return ec;
	}

	virtual void CopyFrom(const MBSMultipleSensor& s)
	{
		Sensor::CopyFrom(s);
		sensors = s.sensors;
		multtype = s.multtype;
		offset = s.offset;
		factor = s.factor;
	}

	virtual mystr GetTypeName()
	{
		return "Multiple sensor";
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer //$ DR 2012-10 return value changed from bool to int

	virtual void AddSensor(int sens_MBS_ID) { sensors.Add(sens_MBS_ID); }

	virtual double GetCurrentValue(double time);

	virtual void SetFactor(double fact) {factor = fact;}
	virtual double GetFactor() const {return factor;}
	virtual void SetOffset(double off) {offset = off;}
	virtual double GetOffset() const {return offset;}


private:
	TMBSMultSensor multtype;
	TArray<int> sensors;

	double factor;
	double offset;
};


#endif