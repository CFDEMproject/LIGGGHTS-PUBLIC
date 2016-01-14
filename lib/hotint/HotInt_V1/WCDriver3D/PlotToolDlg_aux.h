//#**************************************************************
//#
//# filename:             PlotToolDlg_aux.h         
//#
//# author:               Dorninger Alexander
//#
//# generated:						March 2012
//# description:          Classes for own Colorpaletes
//# remarks:						  * originally in PlotTool, 
//#                       * also needed in other routines
//#
//# This file is part of HotInt.
//# HotInt is free software: you can redistribute it and/or modify it under the terms of 
//# the HOTINT license. See folder 'licenses' for more details.
//#
//# bug reports are welcome!!!
//# WWW:		www.hotint.org
//# email:	bug_reports@hotint.org or support@hotint.org
//#***************************************************************************************
 

#ifndef PLOTTOOLDLGAUX__H
#define PLOTTOOLDLGAUX__H

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ classes for Color Palette
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// assign the color a NAME, Vector3D value, intCOLORREF value
// identity operator to search the list by name
// TODO: interface to EDC / textfile
class MyColorInfo
{
public: // lifecylce
	MyColorInfo() { Reset(); }
	MyColorInfo(MyColorInfo& other) { CopyFrom(other); }
	~MyColorInfo() {}
public: // lifecycle II

	MyColorInfo(mystr& namei, Vector3D colorveci, int colorinti ) {  name = namei; colorvec = colorveci; colorint = colorinti; }
	void CopyFrom(MyColorInfo& other)
	{
		name = other.name;
		colorvec = other.colorvec;
		colorint = other.colorint;
	}
	void Reset() { name = mystr(""); colorvec = Vector3D(0.,0.,0.); colorint = 0x00000000; }
public: // access
	mystr& Name() { return name; }
	Vector3D& Color_Vector() { return colorvec; }
	int& Color_COLREF() { return colorint; }
public: // functionality
	int operator == (const MyColorInfo& other) // equality operator for TArray<T>::Find(const T& t)
	{
		// to find in list the name has to match
		return name.Compare(other.name);
	}

	int ReadFromEDC(ElementDataContainer* edc)
	{
    // read all relevant entries that are in the EDC ( less then 3 also valid )
		int flag_name_defined = edc->Find("name");
		if(flag_name_defined)
		{
			Name() = edc->TreeGetString("name");
		}
		int flag_colvec_defined = edc->Find("vec");
		if(flag_colvec_defined)
		{
		  edc->TreeGetVector3D("vec",colorvec.X(),colorvec.Y(),colorvec.Z());
		}
		int flag_colint_defined = edc->Find("rgb");
		if(flag_colint_defined)
		{
			Color_COLREF() = edc->TreeGetInt("rgb");
		}

		if(!flag_name_defined && !flag_colvec_defined && !flag_colint_defined)
			return 0; // error: nothing defined
		if(flag_colvec_defined && !flag_colint_defined)
		{
			int red = (int) (Color_Vector().X()*256. +0.5);
			int green = (int) (Color_Vector().Y()*256. +0.5);
			int blue = (int) (Color_Vector().Z()*256. +0.5);

		  Color_COLREF() = red + 256*green + 256*256*blue;
		}
		if(flag_colint_defined && !flag_colvec_defined)
		{
      int red = (Color_COLREF() & 0x0000FF) >> 0;
			int green = (Color_COLREF() & 0x00FF00) >> 8;
			int blue = (Color_COLREF() & 0xFF0000) >> 16;

			double r = red / 255.;
			double g = green / 255.;
			double b = blue / 255.;
			
			Color_Vector() = Vector3D(r,g,b);
		}
		return 1;
	}

	int WriteToEDC(ElementDataContainer* edc)
	{
		ElementData ed;
		ed.SetText(Name(),"name"); edc->Add(ed);
		ed.SetVector3D(Color_Vector().X(),Color_Vector().Y(),Color_Vector().Z(),"vec"); edc->Add(ed);
		ed.SetInt(Color_COLREF(),"rgb"); edc->Add(ed);
		return edc->Length();
	}


public:
	mystr name;                        // identifier
	Vector3D colorvec;                 // color as Vector3D ( as used for Mesh-, Geom-Objects
	int colorint;                      // color as COLORREF ( as used for PlotTool )
};

// Palette - changable, searchable, default initializaiton
// TODO: interface to EDC / textfile
class MyColorList
{
public: // lifecycle
	MyColorList() { Reset(); AddDefaultColors(); }
	MyColorList(MyColorList& other) { CopyFrom(other); }
	~MyColorList() {}
public: // lifecycle II
	void CopyFrom(MyColorList& other)
	{
		Reset();
		for(int i=1; i<=other.N(); i++)
		{
			Add(other.Item(i));
		}
	}
	void Reset() { the_colors.Flush(); }

public: // access
	MyColorInfo& Item(int i) { return the_colors((i-1)%N()+1); }   // pick item by number - include modulo function
	int GetColRef(int i) { return Item(i).Color_COLREF(); }
	int GetColRef(mystr& name) { return Item(name).Color_COLREF(); }
	Vector3D GetCol3D(int i) { return Item(i).Color_Vector(); }
	Vector3D GetCol3D(mystr& name) { return Item(name).Color_Vector(); }
	mystr& Name(int i) { return Item(i).Name(); }
public: // search array
	MyColorInfo& Item(mystr& name)                           // pick by name
	{ 
		int nr = the_colors.Find( MyColorInfo(name, Vector3D(0), 0) );
		if (nr) return the_colors( nr ); 
		else return the_colors(1);
	}  
	int Find(int colorrefi)                // find by colorref - returns listindex
	{
		for(int i=1; i<N(); i++)
		{
			if( the_colors(i).Color_COLREF() == colorrefi ) 
				return i;
		}
		return 0;
	}
	int Find(Vector3D& colorveci)
	{
		for(int i=1; i<N(); i++)
		{
			if( the_colors(i).Color_Vector() == colorveci ) 
				return i;
		}
		return 0;
	}

public: // array manipulation
	int Add(MyColorInfo& itemi) { return the_colors.Add(itemi); }
	int N() { return the_colors.Length(); }
	int Remove(int i) { the_colors.Erase(i); return N(); }
	int Remove(MyColorInfo& itemi) { the_colors.Erase(the_colors.Find(itemi)); return N();}

public: // initialization - tables and import functions
	void AddDefaultColors(int cutoff = -1)
	{
		Add( MyColorInfo( mystr("black"),	  Vector3D(0.,0.,0.), 0x00000000 ));

		Add( MyColorInfo( mystr("red"),			Vector3D(1.,0.,0.), 0x000000FF ));
		Add( MyColorInfo( mystr("green"),		Vector3D(0.,1.,0.), 0x0000FF00 ));
		Add( MyColorInfo( mystr("blue"),		Vector3D(0.,0.,1.), 0x00FF0000 ));

		Add( MyColorInfo( mystr("magenta"),	Vector3D(1.,0.,1.), 0x00FF00FF ));
		Add( MyColorInfo( mystr("yellow"),	Vector3D(1.,1.,0.), 0x0000FFFF ));
		Add( MyColorInfo( mystr("cyan"),		Vector3D(0.,1.,1.), 0x00FFFF00 ));
		// ADD MORE DEFAULT COLORS HERE...

		if (cutoff > 0)
		{
			while(cutoff < N())
			{			
				Remove(N());
			}
		}
	}

	int ReadFromEDC(ElementDataContainer* edc)
	{
		// disable popup warnings
		int edc_warn_lvl = edc->GetEDCWarningLevel();
		edc->SetEDCWarningLevel(0);

		Reset();
    // process sub-EDCs
		for(int i=1; i<=edc->Length(); i++)
		{
			ElementData& ed = edc->Get(i);
			if(ed.IsEDC())
			{
				mystr subEDCname(ed.GetDataName());
				ElementDataContainer* subEDC = ed.GetEDC();
				if(subEDC)
				{
					MyColorInfo newcolor;
					newcolor.ReadFromEDC(subEDC);      // each color sub-EDC reads itself
					if(newcolor.Name() == mystr(""))
					{
					  newcolor.Name() = subEDCname;
					}
					Add(newcolor);
				}
			}
		}
		edc->SetEDCWarningLevel(edc_warn_lvl);
		return N();
	}
	
	int WriteToEDC(ElementDataContainer* edc)
	{
		for(int i=1; i<=N(); i++)
		{
			ElementDataContainer* newcolor = (ElementDataContainer*) new ElementDataContainer;
			Item(i).WriteToEDC(newcolor);
			ElementData ed;
			ed.SetEDC(newcolor,Item(i).Name());
			edc->Add(ed);
		}
		return edc->Length();
	}

private:
	TArrayDynamic<MyColorInfo> the_colors;
};

// TODO: interface to EDC / textfile
// possible
//root
//{
//	ColorPalette
//	{
//    Black
//		{
//			rgb = 00000000
//			vec = [0.,0.,0.]
//		}
//	}
//}
//
//root
//{
//	ColorPalette
//	{
//    Color1
//		{
//			name = "black"
//			rgb = 00000000
//			vec = [0.,0.,0.]
//		}
//	}
//}
// NAME and either vec or RGB must be set


#endif