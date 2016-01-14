//#***************************************************************************************
//# filename:     elementdataaccess.h
//#
//# author:				Yury Vetyukov, Johannes Gerstmayr
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


#pragma once

class MathFunction;

void SetElemDataVector2D(ElementDataContainer& edc, const Vector2D& v, const char* name, const char* tooltiptext="");
void SetElemDataVector3D(ElementDataContainer& edc, const Vector3D& v, const char* name, const char* tooltiptext="");
void SetElemDataIVector(ElementDataContainer& edc, const IVector& v, const char* name, const char* tooltiptext="");
void SetElemDataVector(ElementDataContainer& edc, const Vector& v, const char* name, const char* tooltiptext="");
void SetElemDataMatrix(ElementDataContainer& edc, const Matrix& v, const char* name, const char* tooltiptext="");
void SetElemDataMathFunc(ElementDataContainer& edc, MathFunction& mathfunc, const mystr& funcname);

//JG2012-01: added new access functions based on ed, not edc:
void SetElemDataVector2D(ElementData& ed, const Vector2D& v, const char* name);
void SetElemDataVector3D(ElementData& ed, const Vector3D& v, const char* name);
void SetElemDataIVector(ElementData& ed, const IVector& v, const char* name);
void SetElemDataVector(ElementData& ed, const Vector& v, const char* name);
void SetElemDataMatrix(ElementData& ed, const Matrix& v, const char* name);

//improved access functions for elementdata
//flag&1 --> needed value
int GetElemDataBool(MBS* mbs, const ElementDataContainer& edc, const char* name, int& v, int flag = 1);
int GetElemDataInt(MBS* mbs, const ElementDataContainer& edc, const char* name, int& v, int flag = 1);
int GetElemDataDouble(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v, int flag = 1);
int GetElemDataText(MBS* mbs, const ElementDataContainer& edc, const char* name, mystr& str, int flag = 1);
int GetElemDataVector2D(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector2D& v, int flag = 1);
int GetElemDataVector3D(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector3D& v, int flag = 1);
int GetElemDataVector2D(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v1, double& v2, int flag = 1);
int GetElemDataVector3D(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v1, double& v2, double& v3, int flag = 1);
int GetElemDataVector(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector& v, int flag = 1);
int GetElemDataIVector(MBS* mbs, const ElementDataContainer& edc, const char* name, IVector& v, int flag = 1);
int GetElemDataMatrix(MBS* mbs, const ElementDataContainer& edc, const char* name, Matrix& v, int flag = 1);
int GetElemDataMathFunc(MBS* mbs, const ElementDataContainer& edc, const mystr& funcname, MathFunction& mathfunc, int flag = 1); 

//find element, assume that it is Vector3D; if not exists use default value
Vector3D EDCTreeGetVector3D(ElementDataContainer& edc, const char* name, Vector3D default_val = Vector3D(0.));
//find element, assume that it is Vector; if not exists use default value
Vector EDCTreeGetVector(ElementDataContainer& edc, const char* name, Vector default_val = Vector(1));

mystr GetRotUnitStr(int type); //0=rad, 1=degree