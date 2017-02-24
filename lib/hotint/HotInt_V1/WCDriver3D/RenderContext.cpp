//#**************************************************************
//# filename:             RenderContext.cpp
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
 

#include "stdafx.h"

#include "RenderContext.h"

int vertexcnt = 0;
int gltrigcnt = 0;
int glquadcnt = 0;
int gllinecnt = 0;

void RenderContext_::ChooseColor(float r, float g, float b)
{
	glColor3f(r,g,b);
}

void RenderContext_::ChooseLineThickness(float thickness)
{
	glLineWidth(thickness);
}

void RenderContext_::ChoosePointSize(float size)
{
	glPointSize(size);
}

#ifdef _DEBUG
void RenderContext_::Begin()
{
	if(bBegin)
		ASSERT(AfxMessageBox("glBegin() called twice without glEnd().\nContinue?",MB_YESNO | MB_ICONSTOP) == IDYES);
	bBegin = true;
}
#endif // _DEBUG

void RenderContext_::glBeginPoints()
{
	Begin();
	glBegin(GL_POINTS);
}

void RenderContext_::glBeginLines()
{
	gllinecnt++;
	Begin();
	glBegin(GL_LINES);
}

void RenderContext_::glBeginLineStrip()
{
	Begin();
	glBegin(GL_LINE_STRIP);
}

void RenderContext_::glBeginLineLoop()
{
	Begin();
	glBegin(GL_LINE_LOOP);
}

void RenderContext_::glBeginTriangles()
{
	gltrigcnt++;
	Begin();
	glBegin(GL_TRIANGLES);
}

void RenderContext_::glBeginTriangleStrip()
{
	Begin();
	glBegin(GL_TRIANGLE_STRIP);
}

void RenderContext_::glBeginTriangleFan()
{
	Begin();
	glBegin(GL_TRIANGLE_FAN);
}

void RenderContext_::glBeginQuads()
{
	glquadcnt++;
	Begin();
	glBegin(GL_QUADS);
}

void RenderContext_::glBeginQuadStrip()
{
	Begin();
	glBegin(GL_QUAD_STRIP);
}

void RenderContext_::glBeginPolygon()
{
	Begin();
	glBegin(GL_POLYGON);
}

void RenderContext_::glEnd()
{
	bBegin = false;
	::glEnd();
}

#ifdef _DEBUG
void RenderContext_::CheckCoordinates(float x,float y,float z)
{
	static bool bKeepCheckingCoordinates = true;
	if(!bKeepCheckingCoordinates)
		return;

	if(fabs(x) > MaxSceneCoord || fabs(y) > MaxSceneCoord || fabs(z) > MaxSceneCoord)
	{
		//would only be necessary if graphics problem
		//ASSERT(AfxMessageBox("The scene coordinates are out of range.\nContinue?",MB_YESNO | MB_ICONSTOP) == IDYES);
		bKeepCheckingCoordinates = false;
	}
}
#endif

void RenderContext_::glVertex(float x, float y, float z)
{
#ifdef _DEBUG
	CheckCoordinates(x,y,-z);
#endif

	vertexcnt ++;
	::glVertex3f(x,y,-z);
}

void RenderContext_::glVertex(const double* v)
{
#ifdef _DEBUG
	CheckCoordinates((float)v[0],(float)v[1],-(float)v[2]);
#endif

	vertexcnt ++;
	::glVertex3f((float)v[0],(float)v[1],-(float)v[2]);
}

void RenderContext_::glRasterPos(float x, float y, float z)
{
	vertexcnt ++;
	//::glRasterPos3f(x,y,-z);
	::glVertex3f(x,y,-z);
}

void RenderContext_::glColor3f(float x, float y, float z)
{
	::glColor3f(x,y,z);
}

void RenderContext_::glColor4f(float x, float y, float z, float transp)
{
	::glColor4f(x,y,z,transp);
}

void RenderContext_::glEnableColorMaterial()
{
	glEnable(GL_COLOR_MATERIAL);
}

void RenderContext_::glDisableColorMaterial()
{
	glDisable (GL_COLOR_MATERIAL);
}

void RenderContext_::glVertex(const Vertex & v)
{
	CheckCoordinates(v.x,v.y,-v.z);
	::glVertex3f(v.x,v.y,-v.z);
}

void RenderContext_::glLighting(int bLightingIsOn)
{
	if(bLightingIsOn) // && pWCDI->GetIOption(206))
		glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
}

void RenderContext_::glNormal(float x, float y, float z)
{
	::glNormal3f(x,y,-z);
}

void RenderContext_::glNormal(const Vertex & v)
{
	::glNormal3f(v.x,v.y,-v.z);
}

void RenderContext_::glNormal(const double* v)
{
	::glNormal3f((float)v[0],(float)v[1],-(float)v[2]);
}

//texture:
void RenderContext_::glInitializeTexture(unsigned int* texName, const void* texImage, int imageH, int imageW) //initialize texture
{
	//glClearColor(0.0, 0.0, 0.0, 0.0);
	//glEnable(GL_DEPTH_TEST);

	if (glstoretextures)
	{

		gltexturecnt++;
		if (gltexturecnt >= nmaxtextures) 
		{
			gltexturecnt = 1;

			if (!gltexturecnt_warned)
			{
				AfxMessageBox("Warning: Maximum number of textures reached!\nTexture will be disrupted\nFurther warnings suppressed ... ", MB_OK | MB_ICONSTOP);
				gltexturecnt_warned = 1;
			}
		}

		*texName = texNames[gltexturecnt];
	}
	else
	{
		*texName = texNames[0]; //one texname has been allocated!
	}

	glBindTexture(GL_TEXTURE_2D, *texName);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); //GL_CLAMP
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); //GL_NEAREST
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); //GL_NEAREST

	if (texturemode == 1 || !glstoretextures) //generate textures
	{
		glTexImage2D(	GL_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texImage);
	}
	else
	{
	  glTexImage2D(	GL_PROXY_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texImage);
	}
}

void RenderContext_::glBeginTexture(const unsigned int& texName) //call this function just before begin QUADS
{
  glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); //GL_DECAL does not support lighting; GL_MODULATE modulates with surface color and supports lighting!
	glBindTexture(GL_TEXTURE_2D, texName);
}

void RenderContext_::glEndTexture()                  //call this function just after glEnd() of QUADS
{
	glDisable(GL_TEXTURE_2D);
}

void RenderContext_::glTexCoord2f(float x, float y)  //local coordinates in texture image (0,0) - (1,1)
{
	::glTexCoord2f(x, y);
}

void RenderContext_::glBeginTexture1D(const unsigned int& texName) //call this function just before begin QUADS
{
  glEnable(GL_TEXTURE_1D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); //GL_DECAL does not support lighting; GL_MODULATE modulates with surface color and supports lighting!
	glBindTexture(GL_TEXTURE_1D, FEcolortexNames[texName]);
}

void RenderContext_::glEndTexture1D()                  //call this function just after glEnd() of QUADS
{
	glDisable(GL_TEXTURE_1D);
}

void RenderContext_::glTexCoord1f(float x) 
{
	::glTexCoord1f(x);
}


void RenderContext_::glPushMatrices()
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
}

void RenderContext_::glPopMatrices()
{
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void RenderContext_::DrawNode(float x, float y, float z)
{
	// TODO
}

void RenderContext_::DrawFixedDof(float x, float y, float z, int dof_num)
{
	// TODO
}

void RenderContext_::DrawForce(float x, float y, float z, float dirx, float diry, float dirz)
{
	// TODO
}




/*
void RenderContext_::glInitializeTexture(unsigned int* texName, const void* texImage, int imageH, int imageW) //initialize texture
{
	//glClearColor(0.0, 0.0, 0.0, 0.0);
	//glEnable(GL_DEPTH_TEST);


	glGenTextures(1, texName);
	glBindTexture(GL_TEXTURE_2D, *texName);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); //GL_CLAMP
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	glTexImage2D(	GL_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, texImage);
}

void RenderContext_::glBeginTexture(const unsigned int& texName) //call this function just before begin QUADS
{
  glEnable(GL_TEXTURE_2D);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D, texName);
}

void RenderContext_::glEndTexture(const unsigned int& texName)                  //call this function just after glEnd() of QUADS
{
	glDisable(GL_TEXTURE_2D);
	glDeleteTextures(1, &texName);
}

*/

int ControlWindowContext_::AddDrawComponent(const DrawComponent& elem)
{
	return the_elements.Add(elem);
}

void ControlWindowContext_::Reset()
{
	the_elements.Flush();
}

int ControlWindowContext_::DrawScene()
{
	return 1;
}

