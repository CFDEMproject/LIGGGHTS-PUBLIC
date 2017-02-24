//#**************************************************************
//# filename:             RenderContext.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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
 

#ifndef __RENDER_CONTEXT_H_INCLUDED__
#define __RENDER_CONTEXT_H_INCLUDED__

#include <math.h>
#include "mbs_interface.h"
#include "..\mbs_interface\rendercontext.h"		// full path - to differ between two headers with the same name
#include "..\WorkingModule\WinCompDriverInterface.h"


// abstract class implementation
class RenderContext_ : protected RenderContext
{	
	bool bBegin;		// is set inside an OpenGL primitive
#ifdef _DEBUG
	void Begin();			// performs the check
#else
	inline void Begin() {}
#endif  // _DEBUG

#ifdef _DEBUG
	void CheckCoordinates(float x,float y,float z);			// performs the check
#else
	inline void CheckCoordinates(float x,float y,float z) {}
#endif  // _DEBUG

	void ChooseColor(float r, float g, float b);
	void ChooseLineThickness(float thickness);
	void ChoosePointSize(float size);

	void glBeginPoints();			// GL_POINTS
	void glBeginLines();			// GL_LINES
	void glBeginLineStrip();		// GL_LINE_STRIP
	void glBeginLineLoop();			// GL_LINE_LOOP
	void glBeginTriangles();		// GL_TRIANGLES
	void glBeginTriangleStrip();	// GL_TRIANGLE_STRIP
	void glBeginTriangleFan();		// GL_TRIANGLE_FAN
	void glBeginQuads();			// GL_QUADS
	void glBeginQuadStrip();		// GL_QUAD_STRIP
	void glBeginPolygon();			// GL_POLYGON
	void glEnd();
	void glVertex(float x, float y, float z);
	void glVertex(const Vertex & v);
	void glRasterPos(float x, float y, float z);
	void glLighting(int bLightingIsOn);
	void glNormal(float x, float y, float z);
	void glNormal(const Vertex & v);
	void glVertex(const double* v);
	void glNormal(const double* v);
	void glColor3f(float x, float y, float z);
	void glColor4f(float x, float y, float z, float transp);
	void glEnableColorMaterial();
	void glDisableColorMaterial();
	void glPopMatrices();
	void glPushMatrices();

	//texture:
	void glInitializeTexture(unsigned int* texName, const void* texImage, int imageH, int imageW); //initialize texture
	void glBeginTexture(const unsigned int& texName); //call this function just before begin QUADS
	void glEndTexture();                  //call this function just after glEnd() of QUADS
	void glTexCoord2f(float x, float y);  //local coordinates in texture image (0,0) - (1,1)


	void glBeginTexture1D(const unsigned int& texName); //call this function just before begin QUADS
	void glEndTexture1D();   //call this function just after glEnd() of QUADS
	void glTexCoord1f(float x);  //local coordinates in texture image (0 - 1)
	
	// defined symbols
	void DrawNode(float x, float y, float z);
	void DrawFixedDof(float x, float y, float z, int dof_num);
								// here dof_num changes from 0 to 2: x, y, z
	void DrawForce(float x, float y, float z, float dirx, float diry, float dirz);
								// dir defines the direction and not the magnitude

	float SetMaxSceneCoord(float msc) {MaxSceneCoord = msc;}
	float GetMaxSceneCoord() const {return MaxSceneCoord;}

protected:

	RenderContext_() : bBegin(false) {}

	float DetailLevel;			// defines current scaling of the scene,
								//   so that the element like force or support can be
								//   drawn in constant size

	float MaxSceneCoord;		// the value returned by WCDInterface::GetSceneMaxAbsCoordinate()

	// the flag, whether to use lighting
	//BOOL bLighting;
	int glstoretextures; //faster redraw, but model needs to be rendered twice!
	int gltexturecnt;
	int gltexturecnt_warned;

	int texturemode; //proxy or not ...
	int nmaxtextures;
	unsigned int* texNames;

  static GLubyte * FEcolortexture; //for 1D - color texture for FE-colors
	unsigned int* FEcolortexNames;
	GLubyte* colortexture;

};

	// Interface to the DrawWindow2D for Control Elements
	//!AD: 2012-12-12 [
class ControlWindowContext_ : public ControlWindowContext
{
// DEFINED IN INTERFACE
	// Add an DrawElement to the List, this is done by the MBS-Elements
public:
	ControlWindowContext_()
	{
		the_elements.Flush();
	}
public:
	int AddDrawComponent(const DrawComponent& elem);

	// DrawFunction for elements of the DrawElementsList
public:
	int DrawScene();

	// reset thhe content of the window
public:
	virtual	void Reset();

// OWN
private:
	TArrayDynamic<DrawComponent> the_elements;
public:	
	DrawComponent& GetElement(int i) { return the_elements(i); }
	int NElements() { return the_elements.Length(); }
	//int DrawComponent(int i);

	//TArray<int> capture;
	//int CaptureScene();
	//int CaptureElement(int i);

	//int DrawRect(int i);


};
//!AD: 2012-12-12 ]

#endif  // __RENDER_CONTEXT_H_INCLUDED__