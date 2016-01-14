//#**************************************************************
//#
//# filename:             PlaneSymmetricTensorComponents.h
//#
//# author:               Vetyukov Yury
//#
//# generated:						February 2011
//# description:          Utility class for ANCFThinPlateElement
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
 
#pragma once

/// a symmetric plane tensor on the surface for the given material coordinates over the element
/// can be defined by three components, which are represented by this class.
/// depending on the context, these are either the covariant or contravariant components.
struct PlaneSymmetricTensorComponents
{
	double t11;
	double t12;
	double t22;
	PlaneSymmetricTensorComponents() : t11(0), t12(0), t22(0) {}
	PlaneSymmetricTensorComponents(double t11, double t12, double t22) { SetComponents(t11, t12, t22); }
	void SetComponents(double t11, double t12, double t22);
	// numbering: alphaBeta = 1: t11, alphaBeta = 2: t12, alphaBeta = 3: t22
	void SetComponent(int alphaBeta, double value);
	// algebraic operations
	PlaneSymmetricTensorComponents & operator*=(double k);
	PlaneSymmetricTensorComponents & operator+=(const PlaneSymmetricTensorComponents & T);
	PlaneSymmetricTensorComponents & operator-=(const PlaneSymmetricTensorComponents & T);
	/// convolution (double inner product) with an
  /// isotropic fourth-rank tensor C: t1**C**t2.
  /// the tensor C is identified by the two constants C1, C2;
  /// the isotropic result depends on the actual metric of the surface,
  /// which is identified by the contravariant components of the first
  /// metric tensor A (in the undeformed configuration).
  /// both t1 and t2 are considered to be covariant components.
  /// in the indexed notation we compute
  /// t1_alpha_beta t2_gamma_delta C^alpha^beta^gamma^delta,
  /// where C^alpha^beta^gamma^delta = C1 A^alpha^beta A^gamma^delta + C2 A^beta^gamma A^alpha^delta.
	double Convolute(const PlaneSymmetricTensorComponents & T, double C1, double C2,
		const PlaneSymmetricTensorComponents & Acontravariant);
	/// computes this**C**this according to Convolute()
	double SelfConvolute(double C1, double C2,
		const PlaneSymmetricTensorComponents & Acontravariant) const;
  /// trace of a tensor, which is computed according to the metric
  /// in the undeformed configuration.
  double Trace(const PlaneSymmetricTensorComponents & Acontravariant);
	/// build the deviator from *this
	/// dev *this = *this - tr(*this)/3 * I
	/// due to the specific metric, the contravariant of A is used for the trace operator,
	/// and the covariant of A replaces the identity operator I
	void BuildDeviator(const PlaneSymmetricTensorComponents & Acovariant, const PlaneSymmetricTensorComponents & Acontravariant);
	/// deteriminant of the matrix of components
  double Det() const;
	/// inverse of the matrix of components - i.e. contravariant components from covariant
  PlaneSymmetricTensorComponents Inverse() const;
};
