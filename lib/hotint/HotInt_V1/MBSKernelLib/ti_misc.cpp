//#**************************************************************
//#
//# filename:             ti_misc.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						10.03.2004
//# description:          Auxiliary Functions for timeint (mainly for graphics)
//# remarks:						  belongs to header file timeint.h
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
 
#include "../WorkingModule/stdafx.h"
#include "ioincludes.h"

#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <math.h>

#include "mbs_interface.h"
#include "rendercontext.h"
#include "timeint.h"

extern UserOutputInterface * global_uo;

//#include "../WorkingModule/stdafx.h"
//#include <gl/gl.h>
//#include <gl/glu.h>


int useFEcolortex = 1;


void TimeInt::SetTexStoreResolution(int n) 
{

	int fact = 2*(int)pow(2.,n);

	if (fact > GetTexStoreMaxResolution()) fact = GetTexStoreMaxResolution();

	texstorenn = fact;
	//uo << "settexstoreres=" << texstorenn << "\n";

	if (texImage != NULL) delete[] texImage;
	texImage = new char[4*fact*fact*4]; //safety factor 4, to ensure that there are no problems ...
}



void TimeInt::SetColor(const Vector3D& col)
{
	pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], 1.f - (1.f-(float)GetDOption(206))*(float)transparency_on);
	//pCurrentRC->ChooseColor((float)col[0],(float)col[1],(float)col[2]);
	actcolor = col;
}

void TimeInt::ChooseColor(float R, float G, float B) const
{
	pCurrentRC->glColor4f(R, G, B, 1.f - (1.f-(float)GetDOption(206))*(float)transparency_on);
	//pCurrentRC->ChooseColor(R, G, B);
}

void TimeInt::SetFEColor(double val) const
{
	if (useFEcolortex)
	{
		if (GetIOption(108) && val >= 0 && val <= 1.) //logarithmic
		{
			//val = 1.-pow(1.-val,10.);
			val = pow(val,0.2);
		}

		/*if (GetIOption(116))  //chladni, 0==black
		{
			if (val < 1e-4) val = 0;
			else val = 1;
		}*/

		if (GetIOption(107)) //invert color
		{
			val = 1.-val; 
		}

		/*
		double niso = GetIOption(103);
		if (niso <= 32) //niso>32 == no tiling
		{
			val = ((int)(niso*val)); //for isolines ...
			val = val/niso+1e-12;
		}*/

	  //call before drawing of object: glEnable(GL_TEXTURE_1D); afterwards glDisable(GL_TEXTURE_1D)!!!
		//if (val <= 0) global_uo << "cval=" << val << "\n";
		//if (val >= 1) global_uo << "cval=" << val << "\n";

		float fval = (float)(0.99*val + 0.005);
		if (val < -1e-8) fval = 0;
		else if (val > 1.00000001) fval = 1;

	  pCurrentRC->glTexCoord1f (fval); //values outside are no problem??
	}
	else
	{
		Vector3D c = FEColor(val);
		ChooseColor((float)c.X(), (float)c.Y(), (float)c.Z());
	}
}


int lastgreyoption=0;
Vector3D TimeInt::FEColor(double val) const
{
	if (GetIOption(108) && val >= 0 && val <= 1.) 
	{
		//val = 1.-pow(1.-val,10.);
		val = pow(val,0.2);
	}

	/*if (GetIOption(116))  //chladni
	{
		if (val < 1e-4) val = 0;
		else val = 1;
	}*/

	if (GetIOption(107)) 
	{
		val = 1.-val; 
	}

	static Vector3D col1(0.25,0.25,0.25);
	static Vector3D col2(0.1,0.1,0.9);
	static Vector3D col3(0.1,0.9,0.9);
	static Vector3D col4(0.1,0.9,0.1);
	static Vector3D col5(0.9,0.9,0.1);
	static Vector3D col6(0.9,0.1,0.1);
	static Vector3D col7(0.7,0.7,0.7);
	static Vector3D col8(0.1,0.1,0.1);

	if (GetIOption(105)) //greymode
	{
		col1=Vector3D(0,0,0);
		col2=Vector3D(0.1,0.1,0.1);
		col3=Vector3D(0.3,0.3,0.3);
		col4=Vector3D(0.5,0.5,0.5);
		col5=Vector3D(0.7,0.7,0.7);
		col6=Vector3D(0.9,0.9,0.9);
		col7=Vector3D(1,0.9,0.9);
		/*if (GetIOption(116))  //chladni
		{
			col6=Vector3D(1,1,1);
			col7=Vector3D(1,1,1);
		}*/
		lastgreyoption = 1;
	}
	else if (lastgreyoption)
	{
		col1=Vector3D(0.25,0.25,0.25);
		col2=Vector3D(0.1,0.1,0.9);
		col3=Vector3D(0.1,0.9,0.9);
		col4=Vector3D(0.1,0.9,0.1);
		col5=Vector3D(0.9,0.9,0.1);
		col6=Vector3D(0.9,0.1,0.1);
		col7=Vector3D(0.7,0.7,0.7);
	}


	double niso = GetIOption(103);
	val = ((int)(niso*val)); //for isolines ...
	val = val/niso+1e-12;


	double fact = 0.25;
	if (val < -1e-10) return col1; 
	if (val < fact)
	{
    return val/fact*col3+(1.-val/fact)*col2;
	}
	val -=fact;
	if (val < fact)
	{
    return val/fact*col4+(1.-val/fact)*col3;
	}
	val -=fact;
	if (val < fact)
	{
    return val/fact*col5+(1.-val/fact)*col4;
	}
	val -=fact;
	if (val < fact*1.0000000001)
	{
    return val/fact*col6+(1.-val/fact)*col5;
	}

  return col7;
}


void TimeInt::UpdateFEMinMaxCol(double val)
{
	if (!GetIOption(101))
	{
		 TImincol = Minimum(TImincol,val);
	}
	else
		TImincol = GetDOption(101);
	if (!GetIOption(100))
	{
		TImaxcol = Maximum(TImaxcol,val);
	}
	else
		TImaxcol = GetDOption(100);
}


void TimeInt::DrawLegend(double yoff) const
{
	int w,h;
	pCurrentRC->GetWindowSize(w,h);

	double h2 = GetIOption(103);

	double sy = (32.*h2)/(double)h;//(double)h*2.;//0.4*768./(double)h;
	double pix_yoff = (8.+16+32*h2+32.*yoff)/(double)h;
	double offx =-1.+16./(double)w;
	double offy = 1.-pix_yoff; //y-range is -1 (bottom) ... +1 (top)
	double sx = 64./(double)w;
	double tile = h2;
	double dsy = sy / tile;
	double val, cv1, cv2;

	pCurrentRC->glPushMatrices();
	ChooseColor(1,0,0);

	for (double i=0; i<1-1e-6; i+=1./tile)
	{
		Vector3D p1(offx,   offy+i*sy,0);
		Vector3D p2(offx+sx,offy+i*sy,0);
		Vector3D p3(offx+sx,offy+dsy+i*sy,0);
		Vector3D p4(offx,   offy+dsy+i*sy,0);

		pCurrentRC->glBeginTexture1D(0); //color texture
		pCurrentRC->glBeginQuads();
		Vector3D n;
		Normal3D(p4,p2,p1,n);

		val = i;
		if (GetIOption(108) && val >= 0 && val <= 1.) 
		{
			val = pow(val,5); //for nonlinear color scale
		}
		Vector3D c1 = FEColor(val);
		cv1 = val;
		val = i+1./tile;

		if (GetIOption(108) && val >= 0 && val <= 1.) 
		{
			val = pow(val,5); //for nonlinear color scale
		}
		Vector3D c2 = FEColor(val);
		cv2 = val;
		//cv1 = (cv1+cv2)*0.5;
		//cv2 = cv1;


		SetFEColor(cv2);
		//pCurrentRC->glColor3f((float)c2[0],(float)c2[1],(float)c2[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		SetFEColor(cv2);
		//pCurrentRC->glColor3f((float)c2[0],(float)c2[1],(float)c2[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		SetFEColor(cv1);
		//pCurrentRC->glColor3f((float)c1[0],(float)c1[1],(float)c1[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		SetFEColor(cv1);
		//pCurrentRC->glColor3f((float)c1[0],(float)c1[1],(float)c1[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);

		pCurrentRC->glEnd();
		pCurrentRC->glEndTexture1D(); //color

		ChooseColor(0.2f,0.2f,0.2f);
		pCurrentRC->ChooseLineThickness(3);
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glEnd();

	}
	Vector3D p1(offx,   offy,0);
	Vector3D p2(offx+sx,offy,0);
	Vector3D p3(offx+sx,offy+sy,0);
	Vector3D p4(offx,   offy+sy,0);

	ChooseColor(0.2f,0.2f,0.2f);
	pCurrentRC->ChooseLineThickness(3);
	pCurrentRC->glBeginLines();
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glEnd();
	pCurrentRC->ChooseLineThickness(1);

	pCurrentRC->glPopMatrices();
}

//interpolate values at quad with local coords 0..1
double QuadInt(double x, double y, double v1, double v2, double v3, double v4)
{
	//order of nodes:
	//
	//  node 4 +------+ node 3    y
	//         |      |           ^
	//         |      |           |
	//         |      |           ---> x
	//  node 1 +------+ node 2
	//
	return (1.-y)*(1.-x)*v1+(1.-y)*(x)*v2+(y)*(x)*v3+(y)*(1.-x)*v4;
}

//interpolate values at quad with local coords 0..1
Vector3D QuadInt(double x, double y, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3, const Vector3D& v4)
{
	//order of nodes:
	//
	//  node 4 +------+ node 3    y
	//         |      |           ^
	//         |      |           |
	//         |      |           ---> x
	//  node 1 +------+ node 2
	//
	return ((1.-y)*(1.-x))*v1+((1.-y)*x)*v2+(y*x)*v3+(y*(1.-x))*v4;
}


int TimeInt::GetTexStoreMaxResolution()
{
	return 128;
};


void TimeInt::DrawIsoQuadTex(Vector3D* p, const Vector3D& n, double* v, int res, char* teximage)
{
	if (res != texstorenn) return; //GetTexStoreResolution()

	//currently only for ISO-Chladni figures
	Vector3D c;

	double nn = texstorenn;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	pCurrentRC->glInitializeTexture(&texName, texImage, texstorenn, texstorenn);

	pCurrentRC->glBeginTexture(texName);
	pCurrentRC->glBeginQuads();

	const float eps = 0.5f/(float)nn;
	float size = (1.f - eps);


	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glTexCoord2f(eps,eps); 
	pCurrentRC->glVertex((float)p[3].X(),(float)p[3].Y(),(float)p[3].Z());

	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glTexCoord2f(eps,size); 
	pCurrentRC->glVertex((float)p[2].X(),(float)p[2].Y(),(float)p[2].Z());

	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glTexCoord2f(size,size); 
	pCurrentRC->glVertex((float)p[1].X(),(float)p[1].Y(),(float)p[1].Z());

	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glTexCoord2f(size,eps); 
	pCurrentRC->glVertex((float)p[0].X(),(float)p[0].Y(),(float)p[0].Z());

	pCurrentRC->glEnd();
	pCurrentRC->glEndTexture();
}

int oldnn = 0;

void TimeInt::DrawIsoQuad(Vector3D* p, const Vector3D& n, double* v)
{
	if (1) //with textures
	{
		if (useFEcolortex)
		{
			pCurrentRC->glBeginTexture1D(0); //color

			Vector3D c;
			pCurrentRC->glBeginQuads();

			SetFEColor(v[3]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[3].X(),(float)p[3].Y(),(float)p[3].Z());

			SetFEColor(v[2]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[2].X(),(float)p[2].Y(),(float)p[2].Z());

			SetFEColor(v[1]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[1].X(),(float)p[1].Y(),(float)p[1].Z());

			SetFEColor(v[0]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[0].X(),(float)p[0].Y(),(float)p[0].Z());

			pCurrentRC->glEnd();
		
			pCurrentRC->glEndTexture1D();
		}
		/*
		else if (GetIOption(116)) //chladni ISO-lines around -/+ 0
		{
			Vector3D c;
			//texture:

			int rn = texstorenn;
			double drn = rn;

			oldnn = texstorenn;

			double nn = texstorenn;
			double intpv;  //interpolate colors

			double chladnirad = 0.+GetDOption(100);
			int chladnisteps = 16;
			double chpx[24], chpy[24];

			if (GetIOption(116)) //chladni ISO-lines around -/+ 0
			{
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//draw chladni figure: iso-line around zero-line
				for (int i=0; i < chladnisteps; i++)
				{
					chpx[i] = chladnirad * sin((double)i/(double)chladnisteps*2.*MY_PI);
					chpy[i] = chladnirad * cos((double)i/(double)chladnisteps*2.*MY_PI);
				}
			}


			for (int i=0; i<rn; i++)
			{
				for (int j=0; j<rn; j++)
				{
					double gx1 = (double)(i)/drn;
					double gy1 = (double)(drn-j-1)/drn;

					if (!GetIOption(116)) //chladni ISO-lines around -/+ 0
					{
						intpv = QuadInt(gx1,gy1,v[0],v[1],v[2],v[3]);
						//c = 128.+127.*sin(gx1*gy1);
						c = 255.*FEColor(intpv);
					} 
					else
					{
						//draw chladni figure: iso-line around zero-line, radius=chladnirad
						double epss = 1e-8;
						double val0 = QuadInt(gx1,gy1,v[0],v[1],v[2],v[3]);
						double val;
						int flag = Sgn(val0);
						c = 255;
						int cnt = 0;
						for (int k=0; k < chladnisteps; k++)
						{
							val = QuadInt(gx1+chpx[k],gy1+chpy[k],v[0],v[1],v[2],v[3]);
							if (flag != Sgn(val) && (fabs(val-val0) > epss))
							{
								cnt++;
								//if (cnt > chladnisteps/3) {c = 0;break;}
							}
						}

						c = 255*(1.f-64.f*((float)cnt - (int)GetDOption(101))/(float)chladnisteps);
						if (c.X() < 0) c=0;
						if (c.X() > 255) c=255;

					}

					texImage[i*texstorenn*4+j*4+0] = (char)c.X();
					texImage[i*texstorenn*4+j*4+1] = (char)c.Y();
					texImage[i*texstorenn*4+j*4+2] = (char)c.Z();
					texImage[i*texstorenn*4+j*4+3] = (char)255;
				}
			}

			pCurrentRC->glInitializeTexture(&texName, texImage, texstorenn, texstorenn);

			pCurrentRC->glBeginTexture(texName);
			pCurrentRC->glBeginQuads();

			const float eps = 0.25f/(float)rn;
			float size = (1.f - eps);
			

			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glTexCoord2f(eps,eps); 
			pCurrentRC->glVertex((float)p[3].X(),(float)p[3].Y(),(float)p[3].Z());

			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glTexCoord2f(eps,size); 
			pCurrentRC->glVertex((float)p[2].X(),(float)p[2].Y(),(float)p[2].Z());

			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glTexCoord2f(size,size); 
			pCurrentRC->glVertex((float)p[1].X(),(float)p[1].Y(),(float)p[1].Z());

			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glTexCoord2f(size,eps); 
			pCurrentRC->glVertex((float)p[0].X(),(float)p[0].Y(),(float)p[0].Z());

			pCurrentRC->glEnd();
			pCurrentRC->glEndTexture();
		}*/ //chladni
		else
		{
			Vector3D c;
			pCurrentRC->glBeginQuads();

			c=FEColor(v[3]);
			pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[3].X(),(float)p[3].Y(),(float)p[3].Z());

			c=FEColor(v[2]);
			pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[2].X(),(float)p[2].Y(),(float)p[2].Z());

			c=FEColor(v[1]);
			pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[1].X(),(float)p[1].Y(),(float)p[1].Z());

			c=FEColor(v[0]);
			pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p[0].X(),(float)p[0].Y(),(float)p[0].Z());

			pCurrentRC->glEnd();
		}
	}
	else
	{
		double nn = 4; //4 is very good
		Vector3D pp[4];
		double vv[4];

		pCurrentRC->glBeginQuads();
		for (double i=1; i<=nn; i++)
		{
			for (double j=1; j<=nn; j++)
			{
				double gx1 = (i-1)/nn;
				double gx2 = (i)*1.001/nn;  
				double gy1 = (j-1)/nn;
				double gy2 = (j)*1.001/nn;  

				vv[0] = QuadInt(gx1,gy1,v[0],v[1],v[2],v[3]);
				vv[1] = QuadInt(gx2,gy1,v[0],v[1],v[2],v[3]);
				vv[2] = QuadInt(gx2,gy2,v[0],v[1],v[2],v[3]);
				vv[3] = QuadInt(gx1,gy2,v[0],v[1],v[2],v[3]);
				
				pp[0] = QuadInt(gx1,gy1,p[0],p[1],p[2],p[3]);
				pp[1] = QuadInt(gx2,gy1,p[0],p[1],p[2],p[3]);
				pp[2] = QuadInt(gx2,gy2,p[0],p[1],p[2],p[3]);
				pp[3] = QuadInt(gx1,gy2,p[0],p[1],p[2],p[3]);

				Vector3D c;
				c=FEColor(vv[3]);
				pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glVertex((float)pp[3].X(),(float)pp[3].Y(),(float)pp[3].Z());

				c=FEColor(vv[2]);
				pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glVertex((float)pp[2].X(),(float)pp[2].Y(),(float)pp[2].Z());

				c=FEColor(vv[1]);
				pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glVertex((float)pp[1].X(),(float)pp[1].Y(),(float)pp[1].Z());

				c=FEColor(vv[0]);
				pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glVertex((float)pp[0].X(),(float)pp[0].Y(),(float)pp[0].Z());

			}
		}
		pCurrentRC->glEnd();
	}
}


//Draw Quads from list of points arranged in a field of n1 x n2 points, Draw color according to value v
void TimeInt::DrawColorQuads(const TArray<Vector3D>& p, const TArray<double>& v, int n1, int n2, 
														 int colormode, int drawlines, int vres)
{
	//uo << "n1=" << n1 << "\n";
	//uo << "n2=" << n2 << "\n";
	//uo << "p.Length=" << p.Length() << "\n";

	//colormode: 0=no color, 1=FE color
	//drawlines: 1=color+field outline, 2=field outline, 3=color+element outline, 4=element outline
	
	//lateron we will store all quads in a list, determine the range and draw it
	if (!GetIOption(100))
	{
		for (int i=1; i<=v.Length(); i++)
			TImaxcol = Maximum(TImaxcol,v(i));
	}
	else
		TImaxcol = GetDOption(100);

	if (!GetIOption(101))
	{
		for (int i=1; i<=v.Length(); i++)
			TImincol = Minimum(TImincol,v(i));
	}
	else
		TImincol = GetDOption(101);

	double diff = TImaxcol-TImincol;

	if (TImincol == 1e10) {TImincol=TImaxcol=0;}
	if (diff == 0) {diff=1;}

	if (IsFlagComputeMinMaxFECol()) return;

	Vector3D n;
	Vector3D c;
	int p1,p2,p3,p4;
	double vcol[4];
	Vector3D piso[4];

	int usenormals = 0;
	int noff = 0;
	if (p.Length() == 2*n1*n2) 
	{
		//uo << "p.Length=" << p.Length() << ", n1*n2=" << n1*n2 << "\n";
		usenormals = 1;
		noff = n1*n2;
	}

	if (drawlines != 2 && drawlines != 4)
	{
		for (int i1=1; i1 < n1; i1++)
		{
			for (int i2=1; i2 < n2; i2++)
			{
				p4=(i2-1)*n1+i1;
				p3=(i2-1)*n1+i1+1;
				p2=(i2)*n1+i1+1;
				p1=(i2)*n1+i1;

				//if (!usenormals) Normal3D(p(p4),p(p2),p(p1),n);
				if (!usenormals) Normal3D(p(p4),p(p3),p(p2),n);

				if (colormode)
				{
						vcol[0] = (v(p1)-TImincol)/diff;
						vcol[1] = (v(p2)-TImincol)/diff;
						vcol[2] = (v(p3)-TImincol)/diff;
						vcol[3] = (v(p4)-TImincol)/diff;

						piso[0] = p(p1);
						piso[1] = p(p2);
						piso[2] = p(p3);
						piso[3] = p(p4);

						DrawIsoQuad(piso, n, vcol);
					/*
					//chladni
					if (GetIOption(116))  //chladni
					{
					if (vres == 1)
					{
					vcol[0] = v(p1);
					vcol[1] = v(p2);
					vcol[2] = v(p3);
					vcol[3] = v(p4);

					piso[0] = p(p1);
					piso[1] = p(p2);
					piso[2] = p(p3);
					piso[3] = p(p4);

					DrawIsoQuad(piso, n, vcol);
					}
					else
					{

					}
					}*/
				}
				else
				{
					//uo << "n1=" << p(p4+noff) << "\n";
					pCurrentRC->glBeginQuads();

					c=FEColor((v(p4)-TImincol)/diff);
					if (colormode) pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
					if (!usenormals)
						pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					else
						pCurrentRC->glNormal((float)p(p4+noff).X(),(float)p(p4+noff).Y(),(float)p(p4+noff).Z());
					pCurrentRC->glVertex((float)p(p4).X(),(float)p(p4).Y(),(float)p(p4).Z());

					c=FEColor((v(p3)-TImincol)/diff);
					if (colormode) pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
					if (!usenormals)
						pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					else
						pCurrentRC->glNormal((float)p(p3+noff).X(),(float)p(p3+noff).Y(),(float)p(p3+noff).Z());
					pCurrentRC->glVertex((float)p(p3).X(),(float)p(p3).Y(),(float)p(p3).Z());

					c=FEColor((v(p2)-TImincol)/diff);
					if (colormode) pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
					if (!usenormals)
						pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					else
						pCurrentRC->glNormal((float)p(p2+noff).X(),(float)p(p2+noff).Y(),(float)p(p2+noff).Z());
					pCurrentRC->glVertex((float)p(p2).X(),(float)p(p2).Y(),(float)p(p2).Z());

					c=FEColor((v(p1)-TImincol)/diff);
					if (colormode) pCurrentRC->glColor3f((float)c[0],(float)c[1],(float)c[2]);
					if (!usenormals)
						pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					else
						pCurrentRC->glNormal((float)p(p1+noff).X(),(float)p(p1+noff).Y(),(float)p(p1+noff).Z());
					pCurrentRC->glVertex((float)p(p1).X(),(float)p(p1).Y(),(float)p(p1).Z());

					pCurrentRC->glEnd();
				}
			}
		}
	}
	int advancedlinemode = 1; //draw nices lines (for presenation), but no redraw during zooming
	if (drawlines == 1 || drawlines == 2)
	{
		//if (drawlines != 2)
		if (advancedlinemode) pCurrentRC->BeginLinesOnSolid();
		int p1,p2;
		double thickness = GetDOption(102);
		for (int i1=1; i1 < n1; i1++)
		{
			p1=i1;
			p2=i1+1;
			MyDrawLine(p(p1),p(p2),thickness);
			p1=(n2-1)*n1+i1+1;
			p2=(n2-1)*n1+i1;
			MyDrawLine(p(p1),p(p2),thickness);
		}
		for (int i2=1; i2 < n2; i2++)
		{
			p1=(i2-1)*n1+1;
			p2=(i2-1)*n1+n1+1;
			MyDrawLine(p(p1),p(p2),thickness);
			p1=(i2)*n1;
			p2=(i2)*n1+n1;
			MyDrawLine(p(p1),p(p2),thickness);
		}
		//if (drawlines != 2) 
		if (advancedlinemode) pCurrentRC->EndLinesOnSolid();
	}
	if (drawlines == 3 || drawlines == 4)
	{
		if (advancedlinemode) pCurrentRC->BeginLinesOnSolid();
		double thickness = GetDOption(102);
		int p1,p2,p3,p4;
		for (int i1=1; i1 < n1; i1++)
		{
			for (int i2=1; i2 < n2; i2++)
			{
				p4=(i2-1)*n1+i1;
				p3=(i2-1)*n1+i1+1;
				p2=(i2)*n1+i1+1;
				p1=(i2)*n1+i1;
				MyDrawLine(p(p1),p(p2),thickness);
				MyDrawLine(p(p2),p(p3),thickness);
				MyDrawLine(p(p3),p(p4),thickness);
				MyDrawLine(p(p4),p(p1),thickness);
			}
		}
		if (advancedlinemode) pCurrentRC->EndLinesOnSolid();
	}
}


void TimeInt::DrawQuad(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4) const
{
	pCurrentRC->glBeginQuads();
	Vector3D n;
	Normal3D(p4,p2,p1,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glEnd();
}

void TimeInt::DrawTrig(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3) const
{
	pCurrentRC->glBeginTriangles();
	Vector3D n;
	Normal3D(p3,p2,p1,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glEnd();
}

void TimeInt::DrawPolygon(const TArray<Vector3D>& p, int drawlines, double linewidth) const
{
	Vector3D n(0.,0.,0.);
	if (p.Length() >= 2)
	{
		Normal3D(p(2),p.Last(),p(1),n);
	}

	pCurrentRC->glBeginPolygon();
	for (int i=1; i <= p.Length(); i++)
	{
		pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
		pCurrentRC->glVertex((float)p(i).X(),(float)p(i).Y(),(float)p(i).Z());
	}
	pCurrentRC->glEnd();

	if (drawlines)
	{
		DrawPolygonOutline(p,linewidth);
	}
}

void TimeInt::DrawPolygonOutline(const TArray<Vector3D>& p, double linewidth) const
{
	Vector3D n(0.,0.,0.);
	if (p.Length() >= 2)
	{
		Normal3D(p(2),p.Last(),p(1),n);
	}
	ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
	pCurrentRC->ChooseLineThickness((float)linewidth);
	pCurrentRC->glBeginLineStrip();

	for (int i=1; i <= p.Length(); i++)
	{
		pCurrentRC->glVertex((float)p(i).X(),(float)p(i).Y(),(float)p(i).Z());
	}
	pCurrentRC->glVertex((float)p(1).X(),(float)p(1).Y(),(float)p(1).Z());

	pCurrentRC->glEnd();
	ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
}

void TimeInt::DrawHex(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4,
		const Vector3D& p5,const Vector3D& p6,const Vector3D& p7,const Vector3D& p8, int drawouterfaces) const
{
	//#outer surface right: 3-4-8-7
	//#outer surface left:  1-2-6-5
	//#bottom: 1-3-4-2
	//#top: 5-7-8-6
	float tn = 0.5;
	if (drawlines == 2) tn = 1.;
	else if (drawlines == 3) tn = 2.;

	pCurrentRC->glBeginQuads();
	Vector3D n;
	Normal3D(p1,p2,p4,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glEnd();

	if (drawlines)
	{
		ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
		pCurrentRC->ChooseLineThickness(tn);
		pCurrentRC->glBeginLineStrip();
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glEnd();
		ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	}

	pCurrentRC->glBeginQuads();
	Normal3D(p5,p8,p6,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p8[0],(float)p8[1],(float)p8[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p6[0],(float)p6[1],(float)p6[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p5[0],(float)p5[1],(float)p5[2]);
	pCurrentRC->glEnd();

	if (drawlines)
	{
		ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
		pCurrentRC->ChooseLineThickness(tn);
		pCurrentRC->glBeginLineStrip();
		pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
		pCurrentRC->glVertex((float)p8[0],(float)p8[1],(float)p8[2]);
		pCurrentRC->glVertex((float)p6[0],(float)p6[1],(float)p6[2]);
		pCurrentRC->glVertex((float)p5[0],(float)p5[1],(float)p5[2]);
		pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
		pCurrentRC->glEnd();
		ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	}

	pCurrentRC->glBeginQuads();
	Normal3D(p1,p7,p5,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p5[0],(float)p5[1],(float)p5[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glEnd();

	if (drawlines || drawlinesh)
	{
		ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
		pCurrentRC->ChooseLineThickness(tn);
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
		pCurrentRC->glEnd();
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p5[0],(float)p5[1],(float)p5[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glEnd();
		ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	}

	pCurrentRC->glBeginQuads();
	Normal3D(p2,p6,p8,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p6[0],(float)p6[1],(float)p6[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p8[0],(float)p8[1],(float)p8[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glEnd();

	if (drawlinesh || drawlines)
	{
		ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
		pCurrentRC->ChooseLineThickness(tn);
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p6[0],(float)p6[1],(float)p6[2]);
		pCurrentRC->glEnd();
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p8[0],(float)p8[1],(float)p8[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glEnd();
		ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	}

	if (drawouterfaces)
	{
	pCurrentRC->glBeginQuads();
	Normal3D(p3,p4,p8,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p8[0],(float)p8[1],(float)p8[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p7[0],(float)p7[1],(float)p7[2]);
	pCurrentRC->glEnd();

	pCurrentRC->glBeginQuads();
	Normal3D(p1,p6,p2,n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p6[0],(float)p6[1],(float)p6[2]);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
	pCurrentRC->glVertex((float)p5[0],(float)p5[1],(float)p5[2]);
	pCurrentRC->glEnd();
	}

}

void TimeInt::DrawCube(const Vector3D& p0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3) const//draws cube with reference point and 3 axes
{
	Vector3D p1 = p0;
	Vector3D p2 = p0+v2;
	Vector3D p3 = p0+v1;
	Vector3D p4 = p0+v1+v2;
	Vector3D p5 = v3+p0;
	Vector3D p6 = v3+p0+v2;
	Vector3D p7 = v3+p0+v1;
	Vector3D p8 = v3+p0+v1+v2;

	DrawHex(p1,p2,p3,p4,p5,p6,p7,p8);
}

//draw a zylinder in 3D, with color specified by value v, endpoints p1 and p2, Radius r, number of surfaces (discretized): tile
void TimeInt::DrawColorZyl(double v, const Vector3D& pz1,const Vector3D& pz2, double rz, int tile)
{
	UpdateFEMinMaxCol(v);
	double diff = TImaxcol-TImincol;

	if (TImincol == 1e10) {TImincol=TImaxcol=0;}
	if (diff == 0) {diff=1;}

	if (IsFlagComputeMinMaxFECol()) return;

	pCurrentRC->glBeginTexture1D(0); //color
	SetFEColor((v-TImincol)/diff);
	DrawZyl(pz1,pz2, rz, tile);
	pCurrentRC->glEndTexture1D();
}

//draw a zylinder in 3D, with endpoints p1 and p2, Radius r, number of surfaces (discretized): tile
void TimeInt::DrawZyl(const Vector3D& pz1,const Vector3D& pz2, double rz, int tile) const
{
	float tn = 1.;
	Vector3D n;
	Vector3D p1(pz1.X()+tidraw_offx,
		pz1.Y()+tidraw_offy,
		pz1.Z()+tidraw_offz);
	Vector3D p2(pz2.X()+tidraw_offx,
		pz2.Y()+tidraw_offy,
		pz2.Z()+tidraw_offz);
	double r = rz;//; 

	//const double pi = 3.14159265358979323;
	double i, phi, dphi;
	double dtile = tile;
	if (dtile == 0) {dtile = 1;}

	Vector3D p3=p1+Vector3D(1,1,1);
	Vector3D v21=p2-p1;
	Vector3D v=v21;
	v.Normalize();
	double np1=v*p1;

	if (v.X() != 0) {p3.X() = 1./v.X()*(np1-p3.Y()*v.Y()-p3.Z()*v.Z());}
	else if (v.Y() != 0) {p3.Y() = 1./v.Y()*(np1-p3.X()*v.X()-p3.Z()*v.Z());}
	else if (v.Z() != 0) {p3.Z() = 1./v.Z()*(np1-p3.Y()*v.Y()-p3.X()*v.X());}
	else {return;}

	Vector3D n1 = p3-p1;
	n1.Normalize();
	Vector3D n2 = v.Cross(n1);
	n2.Normalize();

	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginQuads();
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
		pCurrentRC->glEnd();
	}
	//draw triangles after quads, because of transparency!!!
	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginTriangles();
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glEnd();
	}
	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginTriangles();
		pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
		pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
		pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
		pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
		pCurrentRC->glEnd();

		if (drawlines)
		{
			ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
			pCurrentRC->ChooseLineThickness(tn);
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glEnd();
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
			pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
			pCurrentRC->glEnd();
			ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
		}
	}
}

//draw a zylinder in 3D, with endpoints p1 and p2, Radius r, number of surfaces (discretized): tile
void TimeInt::DrawCone(const Vector3D& pz1,const Vector3D& pz2, double rz, int tile, int drawconelines) const
{
	float tn = 1.;
	Vector3D n;
	Vector3D p1(pz1.X()+tidraw_offx,
		pz1.Y()+tidraw_offy,
		pz1.Z()+tidraw_offz);
	Vector3D p2(pz2.X()+tidraw_offx,
		pz2.Y()+tidraw_offy,
		pz2.Z()+tidraw_offz);
	double r = rz;//; 

	//const double pi = 3.14159265358979323;
	double i, phi, dphi;
	double dtile = tile;
	if (dtile == 0) {dtile = 1;}

	Vector3D p3=p1+Vector3D(1,1,1);
	Vector3D v21=p2-p1;
	Vector3D v=v21;
	v.Normalize();
	double np1=v*p1;

	if (v.X() != 0) {p3.X() = 1./v.X()*(np1-p3.Y()*v.Y()-p3.Z()*v.Z());}
	else if (v.Y() != 0) {p3.Y() = 1./v.Y()*(np1-p3.X()*v.X()-p3.Z()*v.Z());}
	else if (v.Z() != 0) {p3.Z() = 1./v.Z()*(np1-p3.Y()*v.Y()-p3.X()*v.X());}
	else {return;}

	Vector3D n1 = p3-p1;
	n1.Normalize();
	Vector3D n2 = v.Cross(n1);
	n2.Normalize();

	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=p2;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=p2;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginQuads();
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
		pCurrentRC->glEnd();
	}
	//draw triangles after quads, because of transparency!!!
	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginTriangles();
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glEnd();

		if (drawconelines)
		{
			ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
			pCurrentRC->ChooseLineThickness(tn);
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glEnd();
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
			pCurrentRC->glEnd();
			ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
		}
	}
}

void TimeInt::DrawColorArrow(double v, const Vector3D& p1, const Vector3D&p2, double linethickness, double headsize, int resolution)
{
	UpdateFEMinMaxCol(v);
	double diff = TImaxcol-TImincol;

	if (TImincol == 1e10) {TImincol=TImaxcol=0;}
	if (diff == 0) {diff=1;}

	if (IsFlagComputeMinMaxFECol()) return;

	pCurrentRC->glBeginTexture1D(0); //color
	SetFEColor((v-TImincol)/diff);
	DrawArrow(p1,p2, linethickness, headsize, resolution);
	pCurrentRC->glEndTexture1D();
}

void TimeInt::MyDrawArrow(const Vector3D& p1, const Vector3D&p2, const Vector3D& col, 
		double linethickness, double headsize, int resolution) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	ChooseColor((float)col.X(), (float)col.Y(), (float)col.Z());
	DrawArrow(p1, p2, linethickness, headsize, resolution);
}

void TimeInt::DrawArrow(const Vector3D& p1, const Vector3D&p2, double linethickness, double headsize, int resolution) const
{
	Vector3D n;
	double len = (p1-p2).Norm();
	if (len < 1e-10) return;

	if (linethickness == -1) linethickness = len * 0.01;
	if (headsize == -1) headsize = len * 0.04;

	double r1 = linethickness;//'line' thickness
	double r2 = headsize;     //head radius
	double lh = headsize*2; //headlength

	double i, phi, dphi;
	double dtile = resolution;
	if (dtile == 0) {dtile = 1;}

	Vector3D v21=p2-p1;

	Vector3D n1,n2;

	v21.SetNormalBasis(n1,n2);
	Vector3D n3 = -1.*v21;
	n3.Normalize();
	Vector3D n4 = -1.*n3;
	Vector3D p2b = p2 + lh*n3;
	Vector3D v21b = p2b-p1;

	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r1*cos(phi)*n1+r1*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21b;
		Vector3D pz4=r1*cos(phi+dphi)*n1+r1*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21b;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		Vector3D pz5 = r2*cos(phi)*n1+r2*sin(phi)*n2+p2b;
		Vector3D pz6 = r2*cos(phi+dphi)*n1+r2*sin(phi+dphi)*n2+p2b;

		pCurrentRC->glBeginQuads();
		//cylinder:
		pCurrentRC->glNormal(nz1.GetVecPtr());
		pCurrentRC->glVertex(pz1.GetVecPtr());
		pCurrentRC->glNormal(nz2.GetVecPtr());
		pCurrentRC->glVertex(pz4.GetVecPtr());
		pCurrentRC->glNormal(nz2.GetVecPtr());
		pCurrentRC->glVertex(pz3.GetVecPtr());
		pCurrentRC->glNormal(nz1.GetVecPtr());
		pCurrentRC->glVertex(pz2.GetVecPtr());

		//ring:
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz2.GetVecPtr());
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz5.GetVecPtr());
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz6.GetVecPtr());
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz3.GetVecPtr());

		pCurrentRC->glEnd();
	}
	//draw triangles after quads, because of transparency!!!
	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r1*cos(phi)*n1+r1*sin(phi)*n2+p1;
		Vector3D pz4=r1*cos(phi+dphi)*n1+r1*sin(phi+dphi)*n2+p1;

		pCurrentRC->glBeginTriangles();
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(p1.GetVecPtr());
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz4.GetVecPtr());
		pCurrentRC->glNormal(n3.GetVecPtr());
		pCurrentRC->glVertex(pz1.GetVecPtr());
		pCurrentRC->glEnd();
	}
	for (i = 0; i < dtile; i++)
	{
		//triangles at head
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r2*cos(phi)*n1+r2*sin(phi)*n2+p2b;
		Vector3D pz4=r2*cos(phi+dphi)*n1+r2*sin(phi+dphi)*n2+p2b;
		Vector3D nz1 = pz1-p2b;
		Vector3D nz2 = pz4-p2b;
		Vector3D nz3 = 0.5*(nz1+nz2);
		nz1.Normalize();
		nz2.Normalize();
		nz3.Normalize();


		pCurrentRC->glBeginTriangles();
		pCurrentRC->glNormal(nz3.GetVecPtr());
		pCurrentRC->glVertex(p2.GetVecPtr());
		pCurrentRC->glNormal(nz1.GetVecPtr());
		pCurrentRC->glVertex(pz1.GetVecPtr());
		pCurrentRC->glNormal(nz2.GetVecPtr());
		pCurrentRC->glVertex(pz4.GetVecPtr());
		pCurrentRC->glEnd();
	}
}


//draw a zylinder in 3D, with endpoints p1 and p2, Radius r, number of surfaces (discretized): tile
//use tangential vectors pz1dir and pz2dir
void TimeInt::DrawZyl(const Vector3D& pz1, const Vector3D& pz2, 
											const Vector3D& pz1dir, const Vector3D& pz2dir, double rz, int leftend, int rightend, int tile) const
{
	float tn = 1.;
	Vector3D n;
	Vector3D p1(pz1.X()+tidraw_offx,
		pz1.Y()+tidraw_offy,
		pz1.Z()+tidraw_offz);
	Vector3D p2(pz2.X()+tidraw_offx,
		pz2.Y()+tidraw_offy,
		pz2.Z()+tidraw_offz);

	Vector3D dir1(pz1dir);
	Vector3D dir2(pz2dir);
	//dir1 = p2-p1;
	//dir2 = p2-p1;
	dir1.Normalize();
	dir2.Normalize();

	double r = rz;//; 

	//const double pi = 3.14159265358979323;
	double i, phi, dphi;
	double dtile = tile;
	if (dtile == 0) {dtile = 1;}

	Vector3D n1a, n1b, n2a, n2b;
	dir1.SetNormalBasis(n1a,n1b);

	n2a = n1a;
	n2b = n1b;
	dir2.GramSchmidt(n2a);
	n2b = dir2.Cross(n2a);

	//global_uo << "Basis1:" << dir1 << ", " << n1a << ", " << n1b << "=" << n1a.Cross(n1b) << "\n";
	//global_uo << "Basis2:" << dir2 << ", " << n2a << ", " << n2b << "=" << n2a.Cross(n2b) << "\n";
	

	for (i = 0; i < dtile; i++)
	{
		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1a+r*sin(phi)*n1b+p1;
		Vector3D pz2=r*cos(phi+dphi)*n1a+r*sin(phi+dphi)*n1b+p1;
		Vector3D pz3=r*cos(phi)*n2a+r*sin(phi)*n2b+p2;
		Vector3D pz4=r*cos(phi+dphi)*n2a+r*sin(phi+dphi)*n2b+p2;

		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz2-p1;
		Vector3D nz3 = pz3-p2;
		Vector3D nz4 = pz4-p2;
		nz1.Normalize();
		nz2.Normalize();
		nz3.Normalize();
		nz4.Normalize();
		
		pCurrentRC->glBeginQuads();
		pCurrentRC->glNormal((float)nz4[0],(float)nz4[1],(float)nz4[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal((float)nz3[0],(float)nz3[1],(float)nz3[2]);
		pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
		pCurrentRC->glEnd();

		if (leftend)
		{
			Vector3D v = -1*dir1;
			pCurrentRC->glBeginTriangles();
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glEnd();
		}

		if (rightend)
		{
			Vector3D v = dir2;
			pCurrentRC->glBeginTriangles();
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
			pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
			pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
			pCurrentRC->glEnd();
		}

		if (drawlines)
		{
			ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
			pCurrentRC->ChooseLineThickness(tn);
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glEnd();
			pCurrentRC->glBeginLines();
			pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
			pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
			pCurrentRC->glEnd();
			ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
		}
	}
}


void TimeInt::DrawColorSphere(double v, const Vector3D& p, double r, int tile, double fill)
{

	UpdateFEMinMaxCol(v);

	double diff = TImaxcol-TImincol;
	if (diff == 0) {diff=1;}
	if (TImincol == 1e10) {TImincol=TImaxcol=0;}

	if (IsFlagComputeMinMaxFECol()) return;

	pCurrentRC->glBeginTexture1D(0);
	SetFEColor((v-TImincol)/diff);   
	DrawSphere(p, r, tile, fill);
	pCurrentRC->glEndTexture1D();
}

TArray<double> sin_acc_tab_draw_sphere(0);
TArray<double> cos_acc_tab_draw_sphere(0);
TArray<Vector3D> normalvector_acc_tab_draw_sphere(0);

//#define drawsphere_fast
void TimeInt::DrawSphere(const Vector3D& p, double r, int tile, double fill) const
{


	Vector3D pk(p.X()+tidraw_offx,
		p.Y()+tidraw_offy,
		p.Z()+tidraw_offz);

	//const double pi = 3.14159265358979323;
	double i, j, phi, dphi;
	if (tile == 0) {tile = 1;}

#ifdef drawsphere_fast
	if (sin_acc_tab_draw_sphere.Length() != tile+1)
	{
		sin_acc_tab_draw_sphere.SetLen(tile+1);
		cos_acc_tab_draw_sphere.SetLen(tile+1);
		//generate sin/cos table
		for (j = 1; j <= tile+1; j++)
		{
			double phik1 = (j-1)*MY_PI/tile;
			sin_acc_tab_draw_sphere(j) = sin(phik1);
			cos_acc_tab_draw_sphere(j) = cos(phik1);
		}
		normalvector_acc_tab_draw_sphere.SetLen((tile+1)*(tile+1)*4);

		for (j = 0; j < tile; j++)
		{
			for (i = 0; i < tile; i++)
			{
				double cosphik1 = cos_acc_tab_draw_sphere(j+1);
				double cosphik2 = cos_acc_tab_draw_sphere(j+2);
				double sinphik1 = sin_acc_tab_draw_sphere(j+1);
				double sinphik2 = sin_acc_tab_draw_sphere(j+2);

				double cosphi1 = cos_acc_tab_draw_sphere(i+1);
				double cosphi2 = cos_acc_tab_draw_sphere(i+2);
				double sinphi1 = sin_acc_tab_draw_sphere(i+1);
				double sinphi2 = sin_acc_tab_draw_sphere(i+2);

				double xpos1 = r*cosphik1;
				double xpos2 = r*cosphik2;
				double ypos1 = r*sinphik1;
				double ypos2 = r*sinphik2;
			Vector3D n1(0.,1.,0.);
			Vector3D n2(0.,0.,1.);
			Vector3D vr1(xpos1,0,0);
			Vector3D vr2(xpos2,0,0);

				Vector3D nz1(xpos1,ypos1*cosphi1,ypos1*sinphi1); //==>PG normal vector weg
				Vector3D nz2=ypos2*cosphi1*n1+ypos2*sinphi1*n2+vr2;
				Vector3D nz4=ypos1*cosphi2*n1+ypos1*sinphi2*n2+vr1;
				Vector3D nz3=ypos2*cosphi2*n1+ypos2*sinphi2*n2+vr2;

				nz1.Normalize();
				nz2.Normalize();
				nz3.Normalize();
				nz4.Normalize();
				normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*1) = nz1;
				normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*2) = nz2;
				normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*3) = nz3;
				normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*4) = nz4;
			}
		}
	}
#endif

	for (j = 0; j < tile*fill; j++)
	{
		for (i = 0; i < tile; i++)
		{

#ifdef drawsphere_fast
			double cosphik1 = cos_acc_tab_draw_sphere(j+1);
			double cosphik2 = cos_acc_tab_draw_sphere(j+2);
			double sinphik1 = sin_acc_tab_draw_sphere(j+1);
			double sinphik2 = sin_acc_tab_draw_sphere(j+2);

			double cosphi1 = cos_acc_tab_draw_sphere(i+1);
			double cosphi2 = cos_acc_tab_draw_sphere(i+2);
			double sinphi1 = sin_acc_tab_draw_sphere(i+1);
			double sinphi2 = sin_acc_tab_draw_sphere(i+2);

			double xpos1 = r*cosphik1;
			double xpos2 = r*cosphik2;
			double ypos1 = r*sinphik1;
			double ypos2 = r*sinphik2;

			Vector3D pz1(xpos1+pk.X(), ypos1*cosphi1+pk.Y(), ypos1*sinphi1+pk.Z()); //==>PG normal vector weg
			Vector3D pz2(xpos2+pk.X(), ypos2*cosphi1+pk.Y(), ypos2*sinphi1+pk.Z());
			Vector3D pz4(xpos1+pk.X(), ypos1*cosphi2+pk.Y(), ypos1*sinphi2+pk.Z());
			Vector3D pz3(xpos2+pk.X(), ypos2*cosphi2+pk.Y(), ypos2*sinphi2+pk.Z());
			Vector3D nz1 = normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*1);
			Vector3D nz2 = normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*2);
			Vector3D nz3 = normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*3);
			Vector3D nz4 = normalvector_acc_tab_draw_sphere(j+1+i*(tile+1)*4);
#else
			phi = i*2.*MY_PI/tile;
			dphi = 2.*MY_PI/tile;
			double phik1 = j*MY_PI/tile;
			double phik2 = (j+1.)*MY_PI/tile;

			double xpos1 = r*cos(phik1);
			double xpos2 = r*cos(phik2);
			double ypos1 = r*sin(phik1);
			double ypos2 = r*sin(phik2);
			Vector3D n1(0.,1.,0.);
			Vector3D n2(0.,0.,1.);
			Vector3D vr1(xpos1,0,0);
			Vector3D vr2(xpos2,0,0);
			Vector3D pz1=ypos1*cos(phi)*n1+ypos1*sin(phi)*n2+pk+vr1;
			Vector3D pz2=ypos2*cos(phi)*n1+ypos2*sin(phi)*n2+pk+vr2;
			Vector3D pz4=ypos1*cos(phi+dphi)*n1+ypos1*sin(phi+dphi)*n2+pk+vr1;
			Vector3D pz3=ypos2*cos(phi+dphi)*n1+ypos2*sin(phi+dphi)*n2+pk+vr2;
			Vector3D nz1 = pz1-pk;
			Vector3D nz2 = pz2-pk;
			Vector3D nz3 = pz3-pk;
			Vector3D nz4 = pz4-pk;
			nz1.Normalize();
			nz2.Normalize();
			nz3.Normalize();
			nz4.Normalize();
#endif


			pCurrentRC->glBeginQuads();
			pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glNormal((float)nz4[0],(float)nz4[1],(float)nz4[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glNormal((float)nz3[0],(float)nz3[1],(float)nz3[2]);
			pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
			pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
			pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
			pCurrentRC->glEnd();

		}
	}
}

//Draw 2-D line in 3D: (2 coordinates per point ignored)
void TimeInt::MyDrawLineH(const Vector3D& p1, const Vector3D& p2, const Vector3D& vy2, 
													double t, double h, int drawouterface) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	Vector3D v = p2-p1;
	double len = sqrt(sqr(v(1))+sqr(v(2)));
	if (len!=0) v*=1./len;
	len = sqrt(sqr(vy2(1))+sqr(vy2(2)));
	if (len!=0) len=1./len;

	Vector3D nvy1(v(2)*h,-v(1)*h,0);
	Vector3D nvy2(vy2(2)*h*len,-vy2(1)*h*len,0);
	Vector3D nvz(0.,0.,1.*t);
	Vector3D v1(p1(1)+tidraw_offx,p1(2)+tidraw_offy,p1(3));
	Vector3D v2(p2(1)+tidraw_offx,p2(2)+tidraw_offy,p2(3));

	DrawHex(v2-nvy2-nvz,v2-nvy2+nvz,v2+nvy2-nvz,v2+nvy2+nvz,
		v1-nvy1-nvz,v1-nvy1+nvz,v1+nvy1-nvz,v1+nvy1+nvz, drawouterface);

}

void TimeInt::MyDrawLine(const Vector3D& p1, const Vector3D& p2, double thickness) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
	pCurrentRC->ChooseLineThickness((float)thickness);
	pCurrentRC->glBeginLines();
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glEnd();
	ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	pCurrentRC->ChooseLineThickness(1);

}

void TimeInt::MyDrawRectangle(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, double thickness, const Vector3D* colline, const Vector3D* colfill) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	if(colfill)
	{
		ChooseColor((float)(*colfill)[0],(float)(*colfill)[1],(float)(*colfill)[2]);
		pCurrentRC->ChooseLineThickness((float)thickness);
		pCurrentRC->glBeginQuads();

		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);

		pCurrentRC->glEnd();
	}	
	if(colline)
	{
		pCurrentRC->ChooseLineThickness(2.0);
		ChooseColor((float)(*colline)[0],(float)(*colline)[1],(float)(*colline)[2]);
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
		pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);

		pCurrentRC->glEnd();
	}
	ChooseColor((float)actcolor[0],(float)actcolor[1],(float)actcolor[2]);
	pCurrentRC->ChooseLineThickness(1);
}

void TimeInt::MyDrawLine(const Vector3D& p1, const Vector3D& p2, double thickness, const Vector3D& col) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	ChooseColor((float)col[0],(float)col[1],(float)col[2]);
	pCurrentRC->ChooseLineThickness((float)thickness);
	pCurrentRC->glBeginLines();
	pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
	pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
	pCurrentRC->glEnd();
	pCurrentRC->ChooseLineThickness(1);
}

void TimeInt::MyDrawCircleXY(const Vector3D p, double r, const Vector3D& col, int res, double thickness) const
{
	if (IsFlagComputeMinMaxFECol()) return;

	double ns = res;

	//draw points must be Vector3D

	double phi1, phi2;
	for (int i=0; i < ns; i++)
	{
		phi1 = (double)i/ns*2.*MY_PI;
		phi2 = (double)(i+1)/ns*2.*MY_PI;
		
		Vector3D p1(sin(phi1)*r,cos(phi1)*r,0.);
		Vector3D p2(sin(phi2)*r,cos(phi2)*r,0.);

		MyDrawLine(p+p1,p+p2,thickness, col);
	}
}


//+++++++++++++++++++++++++++++++++++++++++
void TimeInt::AddDrawComponent_Line(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D p1, Vector2D p2, Vector3D col)
{
	ControlWindowContext::DrawComponent dcomp;
	dcomp.mbs_type = mbs_type;
	dcomp.mbs_elnr = mbs_elnr;
	dcomp.sub_type = sub_type;
	dcomp.sub_elnr = sub_elnr;
	dcomp.geo_type = TGLine;
	
	dcomp.P1() = p1;
	dcomp.P2() = p2;
  dcomp.Color() = col;
  // dont set 2nd Color
	// dont set text
	p2DDrawWindow->AddDrawComponent(dcomp);
}

void TimeInt::AddDrawComponent_Rect(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background)
{
	ControlWindowContext::DrawComponent dcomp;
	dcomp.mbs_type = mbs_type;
	dcomp.mbs_elnr = mbs_elnr;
	dcomp.sub_type = sub_type;
	dcomp.sub_elnr = sub_elnr;
	dcomp.geo_type = TGRectangle;
	dcomp.Center() = center;
	dcomp.Size() = size;
  dcomp.Color() = col_border;
	dcomp.Background() = col_background;
  // dont set text
	p2DDrawWindow->AddDrawComponent(dcomp);
}

void TimeInt::AddDrawComponent_Ellipse(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background)
{
	ControlWindowContext::DrawComponent dcomp;
	dcomp.mbs_type = mbs_type;
	dcomp.mbs_elnr = mbs_elnr;
	dcomp.sub_type = sub_type;
	dcomp.sub_elnr = sub_elnr;
	dcomp.geo_type = TGEllipse;
	dcomp.Center() = center;
	dcomp.Size() = size;
  dcomp.Color() = col_border;
	dcomp.Background() = col_background;
  // dont set text
	p2DDrawWindow->AddDrawComponent(dcomp);
}

void TimeInt::AddDrawComponent_Text(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_text, mystr& text, TTextAllign positioning)
{
	ControlWindowContext::DrawComponent dcomp;
	dcomp.mbs_type = mbs_type;
	dcomp.mbs_elnr = mbs_elnr;
	dcomp.sub_type = sub_type;
	dcomp.sub_elnr = sub_elnr;
	dcomp.geo_type = TGText;
	dcomp.Center() = center;
	dcomp.Size() = size;
  dcomp.Color() = col_text;
  // dont set 2nd Color
  dcomp.Text() = text;
	dcomp.TextAllign() = positioning;
	p2DDrawWindow->AddDrawComponent(dcomp);
}