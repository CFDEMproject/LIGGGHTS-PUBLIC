//#***************************************************************************************
//# filename:     renderContext.h
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

#pragma once

// PAINTING graphics
struct RenderContext
{
	virtual ~RenderContext() {}

	// the color setting affects all the graphics & text in the graphics window;
	// r, g, b change from 0 to 1
	virtual void ChooseColor(float r, float g, float b) = 0;
	virtual void ChooseLineThickness(float thickness) = 0;	// thickness in pixels
	virtual void ChoosePointSize(float size) = 0;
	struct Vertex
	{
		float x, y, z;
		Vertex() {}
		Vertex(float x_, float y_, float z_) : x(x_),y(y_),z(z_) {}
	};

	// these functions execute corresponding OpenGL commands
	virtual void glBeginPoints() = 0;		// GL_POINTS
	virtual void glBeginLines() = 0;		// GL_LINES
	virtual void glBeginLineStrip() = 0;	// GL_LINE_STRIP
	virtual void glBeginLineLoop() = 0;		// GL_LINE_LOOP
	virtual void glBeginTriangles() = 0;	// GL_TRIANGLES
	virtual void glBeginTriangleStrip() = 0;// GL_TRIANGLE_STRIP
	virtual void glBeginTriangleFan() = 0;	// GL_TRIANGLE_FAN
	virtual void glBeginQuads() = 0;		// GL_QUADS
	virtual void glBeginQuadStrip() = 0;	// GL_QUAD_STRIP
	virtual void glBeginPolygon() = 0;		// GL_POLYGON
	virtual void glEnd() = 0;
	virtual void glColor3f(float x, float y, float z) = 0;
	virtual void glColor4f(float x, float y, float z, float transp) = 0;
	virtual void glEnableColorMaterial() = 0;
	virtual void glDisableColorMaterial() = 0;
	// both commands specify a point of the primitive
	// must be MinCoord <= x,y,z <= MaxCoord, where
	//		MinCoord & MaxCoord have been specified in GetSceneMinMaxCoordinates()
	// the whole figure should be centered around the point (0,0,0),
	//		as the rotation is going around this point
	virtual void glVertex(float x, float y, float z) = 0;
	virtual void glRasterPos(float x, float y, float z) = 0;
	virtual void glVertex(const Vertex & v) = 0;

	// in the mode with lighting the normals start playing a role
	// the normals need not to be normalized
	virtual void glNormal(float x, float y, float z) = 0;
	virtual void glNormal(const Vertex & v) = 0;

	virtual void glVertex(const double* v) = 0;
	virtual void glNormal(const double* v) = 0;

	//use different offset for lines on solid, should improve visibility...
	virtual void BeginLinesOnSolid() = 0;
	virtual void EndLinesOnSolid() = 0;


	//functions to draw texture:
	virtual void glInitializeTexture(unsigned int* texName, const void* texImage, int imageH, int imageW) = 0; //initialize texture
	virtual void glBeginTexture(const unsigned int& texName) = 0; //call this function just before begin QUADS
	virtual void glEndTexture() = 0;                  //call this function just after glEnd() of QUADS
	virtual void glTexCoord2f(float x, float y) = 0;  //local coordinates in texture image (0,0) - (1,1)

	virtual void glBeginTexture1D(const unsigned int& texName) = 0; //call this function just before begin QUADS
	virtual void glEndTexture1D() = 0;   //call this function just after glEnd() of QUADS
	virtual void glTexCoord1f(float x) = 0;  //local coordinates in texture image (0 - 1)



	//If you want to draw without transformation, first call glPushMatrices(), after
	//drawing dont forget to call glPopMatrices()
	virtual void glPushMatrices() = 0;
	virtual void glPopMatrices() = 0;

	// defined symbols
	virtual void DrawNode(float x, float y, float z) = 0;
	virtual void DrawFixedDof(float x, float y, float z, int dof_num) = 0;
	// here dof_num changes from 0 to 2: x, y, z
	virtual void DrawForce(float x, float y, float z, float dirx, float diry, float dirz) = 0;
	// dir defines the direction and not the magnitude

	// the text color is not affected by ChooseColor()
	// here the window coordinates are specified, (-1,-1) -- bottom-left, (1,1) -- top-right
	virtual void PrintText2D(float x, float y, const char * text) = 0;
	// here the structured text can be written. nLineNo >=0;
	// nXPos < 0  -- the text is aligned to the left;
	// nXPos == 0 -- the text is aligned to the center;
	// nXPos > 0  -- the text is aligned to the right;
	virtual void PrintTextStruct(int nLineNo, int nXPos, const char * text) = 0;
	// here the 3D scene coordinates are specified,
	//   so that the text position on the screen is not fixed
	virtual void PrintText3D(float x, float y, float z, const char * text) = 0;

	virtual void SetTextColor(float r, float g, float b) = 0;
	virtual void SetBackgroundColor(float r, float g, float b) = 0;
	virtual void SetCenterOffset(float x, float y, float z) = 0;

	//size of drawing window in pixel
	virtual void GetWindowSize(int& width, int& height) = 0;

	// the following options are set by the user in the driving module
	// in order to tell the working module what to paint and how to paint the scene
	virtual int GetDetailLevel() = 0;	// how many line segments represent an arc
	virtual int	ShowLegend() = 0;
};

// Interface to the DrawWindow2D for Control Elements
//!AD: 2012-12-12 [
struct ControlWindowContext
{
	virtual ~ControlWindowContext() {}

	////// identify the type and number of the element on the MBS-side
	////typedef enum { TNoMBSElem = 0, TIOBlock = 1, TSensor = 2 } TMBSElementType;                          

	////// identify the type and number of the element on the Draw-side
	////typedef enum { TNoSubType = 0, 
	////	     TElementFrame = 1, TElementName = 2,
	////			 TInputNode = 3, TOutputNode = 4,
	////			 TConstructionNode = 5,  TConnectionLine = 6,
	////			 TTextObject = 7,
	////			 TSymbol = 8} 	TDrawElementType;   

	////typedef enum { TNoGeoType = 0, TGLine = 1, TGRectangle = 2, TGEllipse = 3, TGText = 4} TGeomElementType;

	////	typedef enum {HLeft=1, HCenter=5, HRight=9,
	////	VTop=10, VCenter=50, VBottom=90} TTextAllign;

	// this struct holds all informaiton for the object to be drawn
  // each MBS-Side element is mapped to one or more Draw-Elements
	//$ YV 2012-12-21: I'd use Vectors here
	struct DrawComponent
	{
// this defines the element in the MBS - COMBINATION is UNIQUE (for a MBS element...)
		TMBSElementType mbs_type;   // type of the Element ( see enum TMBSElementType )
		int mbs_elnr;               // number of the corresponding element in the MBS
// this defines the component  ( - ? COMBINATION is UNIQUE ? NOT IMPLEMENTED THAT WAY YET) 
		TDrawElementType sub_type;  // type of the Element ( see enum TMBSElementType )
		int sub_elnr;
// this defines the drawn shape
		TGeomElementType geo_type;

		Vector2D center,size;				// positions used for rect, ellipse, line
		Vector3D col,col2;					// foreground and background color of element     (border+area / line / text)
		mystr text;                 // text description - ElementName ...              
    TTextAllign textallign;        // identifier for relative textposition ......   
																																							
// accessfunctions
		Vector2D& Center() { return center; }   
		Vector2D& Size() { return size; }
		Vector2D& P1() { return center; }
		Vector2D& P2() { return size; }
		Vector3D& Color() { return col; } 
		Vector3D& Background() { return col2; }
		mystr& Text() { return text; }
		TTextAllign& TextAllign() { return textallign; }
		

		DrawComponent() : mbs_type(TNoMBSElem), mbs_elnr(0), sub_type(TNoSubType), sub_elnr(0), geo_type(TNoGeoType),
			                center(0), size(0), col(0), col2(0), 
											text(0)
		{
			textallign = TTextAllign (HCenter+VCenter);
		}

	 	DrawComponent(const DrawComponent& other)
		{
			mbs_type = other.mbs_type;
			mbs_elnr = other.mbs_elnr;
			sub_type = other.sub_type;
			sub_elnr = other.sub_elnr;
			geo_type = other.geo_type;

			center = other.center;
			size = other.size;
			col = other.col;
			col2 = other.col2;
			text = other.text;
			textallign = other.textallign;
		}
// Type considering Drawing - Geometry Type
		int IsLine() const { if (geo_type == TGLine) return 1; 	else return 0; }
		int IsRectangle() const { if (geo_type == TGRectangle) return 1; else return 0; }
		int IsEllipse() const { if (geo_type == TGEllipse) return 1; else return 0;	}
		int IsText() const { if (geo_type == TGText) return 1; else return 0;	}

// Tyoe considering Selection - Subtype
		int IsSymbol() const { if (sub_type == TSymbol) return 1; else return 0; }
		
// Type considering Selection		
		int IsSelectable() const { if (IsSelectableMBSElement() || IsSelectableConnectionLine() || IsSelectableConstructionNode()) return 1; else return 0; }
		int IsSelectableMBSElement() const { if (sub_type == TElementFrame || sub_type == TElementName || sub_type == TSymbol) return 1; else return 0; }
		int IsSelectableConnectionLine() const  { if (sub_type == TConnectionLine ) return 1; else return 0; }
		int IsSelectableConstructionNode() const  { if (sub_type == TConstructionNode ) return 1; else return 0; }

		double GetXMin() const
		{
			if (IsLine()) { return Minimum(center.X(),size.X()); }
			else { return center.X()-0.5*abs(size.X()); }
		}
		double GetXMax() const
		{
			if (IsLine()) { return Maximum(center.X(),size.X()); }
			else { return center.X()+0.5*abs(size.X()); }
		}

		double GetYMin() const
		{
			if (IsLine()) { return Minimum(center.Y(),size.Y()); }
			else { return center.Y()-0.5*abs(size.Y()); }
		}
		double GetYMax() const
		{
			if (IsLine()) { return Maximum(center.Y(),size.Y()); }
			else { return center.Y()+0.5*abs(size.Y()); }
		}
	};

	// Add an DrawElement to the List, this is done by the MBS-Elements
	virtual int AddDrawComponent(const DrawComponent& elem) = 0;
	
	// DrawFunction for elements of the DrawElementsList
	virtual int DrawScene() = 0;
	
	// reset thhe content of the window
	virtual	void Reset() = 0;

};
