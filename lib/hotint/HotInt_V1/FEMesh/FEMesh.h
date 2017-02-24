//#**************************************************************
//#
//# filename:             FEMesh.h
//#
//# authors:              Gerstmayr Johannes
//#                       Dorninger Alexander
//#
//# generated:						April 2010
//# description:          Finite Element mesh
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
 

#ifndef FEMESH__H
#define FEMESH__H

#define MESH_STD_TOL 1E-10 // define before include aux
#include "femesh_aux.h"
#include "FEMeshInterface.h"
#include "element.h"
#include "material.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ (abstract base) class FEMesh:  derived: FEMesh2d, FEMesh3d
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class FEMesh : public FEMeshInterface
{
	friend class FEMesh_Generator; // saves many access functions in the class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle
public: 
	FEMesh() : surface((TSetType)TSetFaces,mystr("surface"))
	{ 
	} 
	FEMesh(FEMesh& other)
	{
		CopyFrom(other);
	}
	~FEMesh() 
	{ 
//		ReleaseArray_TemplatePtr(elements);
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle II
public:
	FEMesh(MBS * mbsi) 
	{ 
		mbs = mbsi; 
		thegenerator.SetFEMesh_Generator(mbs,this); 
		InitializeNodeSearchTree();
	}
	virtual void CopyFrom(const FEMesh& e)
	{
// concerning MBS
		mbs = e.mbs;
		settings = e.settings;
		mbsmaterialnumbers = e.mbsmaterialnumbers;
		mbsnodenumbers = e.mbsnodenumbers;
		mbselementnumbers = e.mbselementnumbers;

// mesh data
		points = e.points;
		ReleaseArray_TemplatePtr(elements); // base class TArray<T*>
		for (int i = 1; i <= e.elements.Length(); i++)
		{
			AddElement(*e.elements(i)); // make a copy
		}	
		faces = e.faces;
		areas = e.areas;
		surface = e.surface;

// materials
		ReleaseArray_TemplatePtr(material_data); // base class TArray<*>
		for(int i = 1; i <= e.material_data.Length(); i++)
		{
			AddMaterial(*e.material_data(i));
		}

// boundary conditions
		ReleaseArray_TemplatePtr(load_data);
		for(int i = 1; i <= e.load_data.Length(); i++)
		{
			load_data.Add( (e.load_data(i))->GetCopy() );
		}
		ReleaseArray_TemplatePtr(constraint_data);
		for(int i = 1; i <= e.constraint_data.Length(); i++)
		{
			constraint_data.Add( (e.constraint_data(i))->GetCopy() );
		}
// initial conditions
 		init_disp = e.init_disp;
		init_vel = e.init_vel;	

// quick search and mapping
		ComputeNodeBox();
		InitializeNodeSearchTree();
		ComputeNodesToElements();
		currentselection = e.currentselection;
		mesh_register = e.mesh_register;
		selections = e.selections;

// generator class
		thegenerator = e.thegenerator;
	}	
	// NO GETCOPY - FEMESH IS AN ABSTRACT CLASS
	virtual void Reset()
	{
		mbsmaterialnumbers.Flush();
		mbsnodenumbers.Flush();
		mbselementnumbers.Flush();

		points.Flush();
		ReleaseArray_TemplatePtr(elements);
		faces.Flush();
		areas.Flush();
		surface.FlushArrays();

		ReleaseArray_TemplatePtr(material_data);

		ReleaseArray_TemplatePtr(load_data);
		ReleaseArray_TemplatePtr(constraint_data);

		init_disp.Flush();
		init_vel.Flush();

		nodebox.Clear();
		nodestree.ClearItems();
		ReleaseArray_TemplatePtr(nodes_to_elements);
		currentselection.FlushArrays();
		mesh_register.FlushArrays();
		selections.Flush();
	}

private:
	virtual void Abstractor() =0; // this is an abstract base class !

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ access - prototypes - to be overwritten in derived class
public:
	virtual int GetMeshDim() const { return 0; } // returns mesh dimension

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ access and manipulation of arrays

// concerning MBS
public:	
	// $EK 2013-02-18 add non-const GetMBS() ... needed for GetMBS()->UO()
	virtual MBS * GetMBS() {return mbs;} // returns pointer to mbs
	virtual const MBS * GetMBS() const {return mbs;} // returns pointer to mbs
	FEMesh_Settings& Settings() const { return settings; }		// this function is declared as "const" to avoid changes in the client code
	virtual int GetMBSMaterialNumber(int meshmaterialnumber) { return mbsmaterialnumbers(meshmaterialnumber); }
	virtual const int GetMBSMaterialNumber(int meshmaterialnumber) const { return mbsmaterialnumbers(meshmaterialnumber); }
	virtual int GetMBSNodeNumber(int meshnodenumber) { return mbsnodenumbers(meshnodenumber); }
	virtual const int GetMBSNodeNumber(int meshnodenumber) const { return mbsnodenumbers(meshnodenumber); }
	virtual int GetMBSElementNumber(int meshelementnumber) { return mbselementnumbers(meshelementnumber); }
	virtual const int GetMBSElementNumber(int meshelementnumber) const { return mbselementnumbers(meshelementnumber); }

// mesh data - nodes
// NOTE: data contained in derived class
// THINK ABOUT: back to base class =
// todo: get rid of *POINT* functions
public:
	virtual TArray<FEMesh_Node>& GetNodesArray() { return points; }
	virtual int NNodes() const { return points.Length(); }
	virtual int NPoints() const { return points.Length(); } 
	virtual int NN() const { return points.Length(); } 
	
	virtual const FEMesh_Node& GetNode(int i) const { return points(i); }
	virtual FEMesh_Node& GetNode(int i) { return points(i); }

	virtual const Vector3D GetPoint3D(int i) const { return points(i).GetCoords3D(); } // return 3D vector of specified nodes coordinates	
	virtual const Vector2D GetPoint2D(int i) const { return points(i).GetCoords2D(); } // return 2D vector of specified nodes coordinates
	virtual const int GetNodeDomain(int i) const { return points(i).Domain(); } // returns domain number of specified node
	
	virtual int SetPoint(const int i, const Vector3D p) { points(i).SetCoords3D(p); return 0; } // sets coorninates for a specific node
	virtual int SetPoint(const int i, const Vector2D p) { points(i).SetCoords2D(p); return 0; } // sets coorninates for a specific node
	virtual int SetPointsArraySize(const int i, int flag_clear = 0)  	// reset size of points array, optional reset of entries
	{ 
		points.SetLen(i);
		if (flag_clear) points.SetAll(FEMesh_Node());
		return points.Length();
	}

	virtual int AddPoint(const Vector3D p, int domain = 1) { return points.Add(FEMesh_Node(p,domain)); } // add a node to the node list in mesh
	virtual int AddPoint(const Vector2D p, int domain = 1) { return points.Add(FEMesh_Node(p,domain)); } // add a node to the node list in mesh

	virtual int AddNodeCheck(Vector3D pos, int domain = 1, double tol=MESH_STD_TOL); // creates node if not already present 
	virtual int AddNodeCheck(Vector2D pos, int domain = 1, double tol=MESH_STD_TOL); // creates node if not already present
	

// mesh data - elements
public:
	virtual TArray<FEElement*>& Elements() { return elements; }
	virtual int NElements() const { return elements.Length(); } // returns number of elements
	virtual const FEElement* ElementPtr(int i) const { return elements(i); } // returns pointer to specified element
	virtual FEElement* ElementPtr(int i) { return elements(i); } // returns pointer to specified element
	virtual const FEElement& GetElement(int i) const { return *elements(i); } // returns reference to specified element
	virtual FEElement& GetElement(int i) { return *elements(i); } // returns reference to specified element
	virtual int AddElement(const FEElement& fe) {	FEElement* elem = fe.GetCopy();	return elements.Add(elem); } //add copy of element to the element list in mesh (OLD)
	virtual int AddElement(FEElement* fep) { return elements.Add(fep); } //add element (directly) to the elementlist in mesh

// mesh data - faces
public:
	virtual const TArray<FEMesh_Face>& Faces() const { return faces; }
	virtual TArray<FEMesh_Face>& Faces() { return faces; }
	virtual int NFaces() const { return Faces().Length(); }
	virtual const FEMesh_Face& GetFace(int i) const { return Faces()(i); }
	virtual FEMesh_Face& GetFace(int i) { return Faces()(i); }
	virtual int AddFace(const FEMesh_Face& face) { return Faces().Add(face); }
	virtual int GetElementOfFace(int i) { return GetFace(i).GetElement(); }
	virtual int GetElementSideOfFace(int i) { return GetFace(i).GetSide(); }

// mesh data - areas
public:
	virtual const TArrayDynamic<FEMesh_Set>& Areas() const { return areas; }
	virtual TArrayDynamic<FEMesh_Set>& Areas() { return areas; }
	virtual int NAreas() const { return Areas().Length(); } // returns number of areas
	virtual const FEMesh_Set& GetArea(int i) const {return Areas()(i);} // returns reference to specified area
	virtual FEMesh_Set& GetArea(int i) {return Areas()(i);} // returns reference to specified area
	virtual int AddArea(const FEMesh_Set& area) { return Areas().Add(area);	} //add copy of area to the area list
  virtual void AddFaceToArea(int fnr, int anr) { GetArea(anr).Faces().Add(fnr); } //add face number to area by number

// mesh data surface
public:
	virtual FEMesh_Set& Surface() { return surface; }

// materials
public:
	virtual int NMaterials() const {return material_data.Length();} // returns number of materials
	virtual void GetMaterial(int matnum, double& rho, double& Em, double& nu) const; // returns parameters of a specified material
	virtual void GetMaterial(int matnum, double& rho, double& Em, double& nu, Vector3D& col) const; // returns parameters of a specified material
	virtual Material& GetMaterial(int matnum) {return *material_data(matnum);} // returns material
	virtual Material* GetMaterialPtr(int matnum) {return material_data(matnum);} // returns pointer to material 
	virtual int AddMaterial(double rho, double Em, double nu, Vector3D color = Vector3D(1.0,0.0,0.0));
	virtual int AddMaterial(Material& mat);
	virtual void ReplaceMaterial(int i, Material& mat);
	virtual void ReplaceMaterial(int i, Material* p_mat);
	virtual int AddElasticMaterial(double rho, double Em, double nu); // add material to the material list(s), if it does not exist:

// boundary conditions - loads 
public:
	virtual int NLoads() const { return load_data.Length(); } // returns number of loads
	virtual void LoadsArrayFlush() { load_data.Flush();} // reset all loads
	virtual FEMesh_Load& GetLoad(int i) { return *load_data(i); } // returns reference to a specified load
	virtual int AddLoad(FEMesh_Load& loadi) // adds load to the loads array
	{
		FEMesh_Load* load = loadi.GetCopy();
		return load_data.Add(load);
	}

// boundary conditions - constraints
public:
	virtual int NConstraints() const { return constraint_data.Length(); } // returns number of constraints
	virtual void ConstraintsArrayFlush() { constraint_data.Flush(); } // reset all constraints
	virtual FEMesh_Constraint& GetConstraint(int i) { return *constraint_data(i); } // returns reference to a specified constraint
	virtual int AddConstraint(FEMesh_Constraint& constrainti) // adds constraint to the constratins array
	{ 
		FEMesh_Constraint* constraint = constrainti.GetCopy();
		return constraint_data.Add(constraint); 
	}

// initial conditions
	virtual void InitializeInitDisplacements() { init_disp.SetLen(NN()); init_disp.SetAll(Vector3D(0.,0.,0.)); } //set length of array to number of nodes and all initial displacements to zero
	virtual void InitializeInitVelocities() { init_vel.SetLen(NN()); init_vel.SetAll(Vector3D(0.,0.,0.)); } //set length of array to number of nodes and all initial velocities to zero

	virtual Vector3D GetInitialNodalDisplacement3D(int i) { return init_disp(i); } // return 3D vector of specified nodes initial displacement
	virtual Vector2D GetInitialNodalDisplacement2D(int i) { return Vector2D(init_disp(i).X(),init_disp(i).Y()); } // return 2D vector of specified nodes initial displacement
	virtual Vector3D GetInitialNodalVelocity3D(int i) { return init_vel(i); } // return 3D vector of specified nodes initial velocity
	virtual Vector2D GetInitialNodalVelocity2D(int i) { return Vector2D(init_vel(i).X(),init_vel(i).Y()); } // return 2D vector of specified nodes initial velocity

	virtual void SetInitialNodalDisplacement(int i, const Vector3D& disp) { init_disp(i) = disp; } // set initial displacement for specified node
	virtual void SetInitialNodalDisplacement(int i, const Vector2D& disp) { init_disp(i) = Vector3D(disp.X(),disp.Y(),0.0); } // set initial displacement for specified node
	virtual void SetInitialNodalVelocity(int i, const Vector3D& vel) { init_vel(i) = vel; } // set initial velocity for specified node
	virtual void SetInitialNodalVelocity(int i, const Vector2D& vel) { init_vel(i) = Vector3D(vel.X(),vel.Y(),0.0); } // set initial velocity for specified node
	
	virtual void AddInitialTranslationalVelocity(const Vector3D& vel);  //add inital translational velocity 'vel' to all nodes
	virtual void AddInitialTranslationalVelocity(const Vector2D& vel) { AddInitialTranslationalVelocity(vel.MakeV3D()); }  //add inital translational velocity 'vel' to all nodes
	virtual void AddInitialRotationalVelocity(const Vector3D& omega, const Vector3D& point_at_rot_axis); 	//add initial rotational velocity field (w)x(r'-r) to all nodes
	virtual void AddInitialRotationalVelocity(double omegaz, const Vector2D& point_at_rot_axis) { AddInitialRotationalVelocity(Vector3D(0.,0.,omegaz), point_at_rot_axis.MakeV3D());}  	//add initial rotational velocity field (w)x(r'-r) to all nodes


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ quick search and mapping
public:
	virtual void ComputeNodeBox(IVector& subset = IVector()); // computes box that contains all nodes
	virtual Box3D& GetNodeBox() {return nodebox;} // returns nodebox
	
	virtual void InitializeNodeSearchTree(); // initializes searchtree (sets bins, adds nodes) with current nodes in the mesh
  virtual void ResizeNodeSearchTree(Box3D addbox, int addnodes); // resizes the Searchtree (adds new region) - used in GenerateBlock
	virtual int3 GetSearchTreeDivisions(Box3D& box, int nnodes); // compute divisions for the Searchtree, cells should be roughly cubes
	virtual SearchTree& NodesTree() { return nodestree; } // access the searchtree containing all nodes
	
	virtual void ComputeNodesToElements(); // generates a list that holds a list of elements the node is part for each node
	virtual TArray<TArray<int>*>& NodesToElementList(int flag_autorecompute = 1); // access node_to_element_list, automatic recomppute if number of entries differs from number of nodes  
	virtual TArray<int>& GetElementsOfNode(int node) { return *(NodesToElementList()(node)); } // returns a specified node's elementlist
	virtual int GetFirstElementOfNode(int node) { return NodesToElementList()(node)->Get(1); } // returns first entry of that node's elementlist

	virtual FEMesh_Set& CurSel() { return currentselection; } // access to selection buffer
	virtual int& CurSel(int i) { return CurSel()(i); } // access to selection buffer

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ generator class
public:
	virtual const FEMesh_Generator& Generator() const { return thegenerator; }
	virtual FEMesh_Generator& Generator() { return thegenerator; }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ indirect - statistics - these quantities are computed
public: 
	virtual int NBodies() { return mesh_register.NBodies(); }
	virtual int CountBodies(IVector& subset = IVector()); // checks all used bodynumbers - updates mesh_register

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ import functions - read from external files

// general import functions
public: 
	virtual int LoadNodes(const mystr & filename, const mystr & format = mystr("ANSYS") );	// loads Nodes data from various sources
	virtual int LoadElements(const mystr & filename, const mystr & format = mystr("ANSYS"), const mystr & writeto = mystr("elements") ); // loads Elements data from various sources, writes them into chosen Array (which must be defined in FEMesh)
	virtual int LoadMaterials(const mystr & filename);																			// loads Materials data from File - replaces any previous entries in materials_* arrays
	virtual int LoadMaterialsEDC(const mystr & filename);																		// loads Materials data from File, assuming EDC structure there - replaces any previous entries in materials_* arrays //$ DR 2012-10

	virtual int LoadLoads(const mystr & filename, const mystr & format = mystr("FORMAT1"));	// loads Loads data - NOT IMPLEMENTED YET
	virtual int LoadConstraints(const mystr & filename,const mystr & format = mystr("FORMAT1"));// loads Constraints data - NOT IMPLEMENTED YET
	virtual int LoadFaces2Areas(const mystr & filename);																		// loads a mapping that combines several faces to areas 

// source specific calls Ansys
public:
	virtual int LoadAnsys2HotintNODES(const mystr & filename);		// load Nodes data (generated with Ansys2Hotint.mac - BlockTag: [NODES])
	virtual int LoadAnsys2HotintELEMENTS(const mystr & filename); // load Elements data (generated with Ansys2Hotint.mac - BlockTag: [ELEMENTS])
	virtual int LoadAnsys2HotintFACES(const mystr & filename); 		// load Faces data (generated with Ansys2Hotint.mac - BlockTag: [FACES])
	virtual int LoadLoads_FORMAT1(const mystr & filename);				// loads Loads data - Text file is in a specific format - see documentation

// import helper functions: consistent node numbering within the elements
public:
	virtual int MapNodesAnsys2Hotint(int elementdimension, int elementorder, int nodesused, TArray<int>& nodeshotint); // changes node order (ansys) to be node order (hotint)
	virtual int MapNodesHotint2Ansys(int elementdimension, int elementorder, int nodesused, TArray<int>& nodesansys); // changes node order (hotint) to be node order (ansys)


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ forward data to MBS

	virtual void AddMeshToMBS(int bodyindex = -1) // bodyindex is obsolete for mesh
	{
		if(this->nodes_to_elements.Length() != NN())
			ComputeNodesToElements();
		AddMaterialsToMBS();
		AddNodesToMBS(bodyindex);
		AddElementsToMBS();
		AddLoadsToMBS();
		AddConstraintsToMBS();
	}
// TODO: declare protected to find direct calls in model files...
// forward by entity type
public:
	virtual void AddPointsToMBS(int bodyindex); //pass node data from internal array to mbs
	virtual void AddNodesToMBS(int bodyindex) {return AddPointsToMBS(bodyindex);} //pass node data from internal array to mbs
	virtual int AddSingleNodeToMBS(int nodenr, SearchTree& addtree, int forced_nodedim = 0); // adds a single node to the MBS, node dimension is computed from existing elements if not forced
	virtual void AddMaterialsToMBS();	//pass material data from internal array to mbs
	virtual void AddElementsToMBS(); //pass element data from internal array to mbs, !nodes and materials must be already added!
	virtual void AddLoadsToMBS(); //pass load data from internal array to mbs, !elements must be already added!
  virtual void AddConstraintsToMBS(); // pass constraint data from internal array to mbs, !elements must be already added!

// helper functions - elements
public:
	// returns new MBSElement of FEelement "elem_num"
	// called by: FEMesh::AddElementsToMBS & CMSElement::AddFEMesh 
	// uses settings (changed by YV, 03.11.10)
	virtual Element* GetPtrToMBSElement_new(int elem_num) const; 

// helper functions - loads
protected:
	virtual void AddLoadToMBS_TAreaLoad(FEMesh_AreaLoad& load); // Processes TAreaLoad type Load and adds to MBS
	virtual void AddLoadToMBS_TFaceLoad(FEMesh_FaceLoad& load); // Processes TFaceLoad type Load and adds to MBS
	virtual void AddLoadToMBS_TBodyLoad(FEMesh_BodyLoad& load); // Processes TBodyLoad type Load and adds to MBS
	virtual void ComputeNodeWeight(TArray<int>& facelist, TArray<double>& nodes_weight, TArray<int>& elemnrs, TArray<int>& local_node_nr); // computes the relative weight of each node in a TAreaLoad or TFaceLoad 
	virtual void AddNodalLoadToMBS(int nodenumber, Vector3D& loadvector, int elemnr, int local_node_nr, TArray<StepSettings>& loadsteps = TArray<StepSettings>(0)); // adds the nodal load to the MBS 

// helper functions - constraints
protected:
	virtual void AddConstraintToMBS_TAreaConstraint(FEMesh_AreaConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr); // Processes TAreaConstraint type Constraint and fill nodeconstraint list
	virtual void AddConstraintToMBS_TFaceConstraint(FEMesh_FaceConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr); // Processes TFaceConstraint type Constraint and fill nodeconstraint list
	virtual void AddConstraintToMBS_TNodalConstraint(FEMesh_NodeConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr); // Processes TNodalConstraint type Constraint and fill nodeconstraint list
	virtual void AddNodeConstraintToMBS(int nodenumber, int direction, int penaltyflag, Vector3D& stiffness, Vector3D& damping, TArray<StepSettings>& steps = TArray<StepSettings>(0)); // adds the nodal constraint to the MBS
	
	virtual void AddConstraintToMBS_TAreaContact_Filter(FEMesh_AreaContact& constraint, TArray<TArray<int>*>& sphericaljointlist, int nr); // Processes TAreaContact type Constraint and fill spherical joint list
	virtual void AddConstraintToMBS_TAreaContact(FEMesh_AreaContact& constraint, TArray<TArray<int>*>& sphericaljointlist, int nr); // Processes TAreaContact, adds sperical joints to mbs (uses sphericaljointlist to remove redundant entries)


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ export functions - export in specified format

// export surface in STL format - can be read to GeomMesh3D object
public:
	virtual void ExportSurfaceTrigsToSTL(char* file, int binary = 0); 
	virtual int SurfaceAsTrigs(TArray<int3>& surfacetrigs); // create a list of all surfaceelement as trigs ( nodenumbers )
	virtual int SurfaceTrigsAsPointSequenceTArray(TArray<Vector3D>& pointsequence); // creates an TArray list of surface trigs defined by 3 consecutive vector3d
	virtual void SurfaceTrigsAsPointSequence(Vector3D** points, int& len); // creates an array of surface trigs defined by 3 consecutive vector3d
	virtual void SurfaceTrigsAsDoubleSequence(double** doubles, int& len); // creates an array of surface trigs defined by 9 consecutive doubles (3 points x 3 coordinates)

// various export formats
public:
//$ YV 2011-03-17: Export of elements, nodes and material properties to a simple text file
	bool ExportAsText(const char * filename);
//$ AD 2012-06-01: Export nodes and elements
	bool ExportNodesAndElements(const char * filename,const mystr& format=mystr("ANSYS"));
//$ LA 2011-06-09: Export of geometry to ABAQUS input file - INCOMPLETE
	bool ExportToAbaqus(const char * filename);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ mesh opetations - 

// elements
	virtual int FindHighestMaterialNumber(FEMesh_Set subset_elements = FEMesh_Set());

// faces (and areas)
public:
	virtual int MapFacesToElements(FEMesh_Set& subset_faces = FEMesh_Set(), int trybothcycles = 0); 	//finds elementnumber and sidenumber for a subset of faces
	virtual int MapAreasToElements(TArray<int> & areastomap, int trybothcycles = 0); //finds elementnumber and sidenumber for a subset of areas
	virtual int FindElementAndSideForFEFace(int facenumber, int trybothcycles = 0); //finds elementnumber and sidenumber for a FEMesh_Face with only nodes known
	virtual int FacesFromNodes(TArray<int>& facenrs, TArray<int>& nodenrs, TArray<int>& subset_validelements = TArray<int>()); // create faces from a list of nodes (can be restricted to subset of elements)
	virtual int FacesFromNodes(FEMesh_Set& rv_faces, FEMesh_Set& in_nodes, FEMesh_Set& subset_validelements = FEMesh_Set())
	{
		return FacesFromNodes(rv_faces.Faces(), in_nodes.Nodes(), subset_validelements.Elements());
	}

// surface
public:
	virtual void FindSurface(FEMesh_Set& subset_of_any_type = FEMesh_Set()); 	//creates all faces of the surfce in list FACES, saves according numbers in set SURFACE, optional filter NOT WORKING YET
	virtual void GetNeighbourElements(int elem, TArray<int>& neighbours); //get a list of elements which are neighbours of this element 
	virtual void GetNeighbourElements_Side(int elem, TArray<int> &neighbours); //get a list of elements which are side-neighbours of this element (subset of node neighbors)
	virtual int FaceInNeighbour(TIntX& face, const TArray<int>& neighbours); //returns neighbour element number for the face or 0	

// mirror
//TODO: mirror only a subset
public:
	virtual void MirrorMeshAtPlane(Vector3D nplane, double cplane, double tol = MESH_STD_TOL, int flag_recalc_surface = 1); //adds new nodes, elements, faces, areas and new loads/constraints to internal arrays 	//calls FindSurface() unless flag is set to 0
	virtual void MirrorMeshAtPlane(Vector2D nplane, double cplane, double tol = MESH_STD_TOL, int flag_recalc_surface = 1); //adds new nodes, elements, faces, areas and new loads/constraints to internal arrays 	//calls FindSurface() unless flag is set to 0
protected:
	void DuplicateNodesAtPlane(IVector& mirrornode, Vector3D nplane, double cplane, double tol = MESH_STD_TOL); 	// duplicates Nodes for mirrormesh
	void DuplicateExistingElements(IVector& mirrorelement, const IVector& mirrornodes);	// duplicates Elements for mirrormesh - assume Nodes already copied
	void DuplicateExistingFaces(IVector& mirrorfaces, const IVector& mirrornodes, const IVector& mirrorelements);	// duplicates Faces for mirrormesh - assume Nodes&Elements already duplicated - only faces that have a mapping to a element/side
	void DuplicateExistingAreas(IVector& mirrorareas, const IVector&mirrorfaces);	// duplicate Areas for mirrormesh - assume faces already duplicated
	void DuplicateExistingLoads(const IVector& mirrornodes, const IVector& mirrorelements, const IVector& mirrorfaces, const IVector& mirrorareas, const Vector3D& nplane);	// duplicates existing Loads - assume nodes, elements, faces and areas duplicated
	void DuplicateExistingConstraints(const IVector& mirrornodes, const IVector& mirrorelements, const IVector& mirrorfaces, const IVector& mirrorareas, const Vector3D& nplane);	// duplicates existing Constraints  - assume nodes, elements, faces and areas duplicated

// refine
public:
	virtual void RefineMesh2(IVector& subset = IVector()); // ATTENTION: works with hexahedrals only so far // refine the mesh with a factor 2 in all directions 
protected:
	virtual void RefinedNodes(IVector& points); // ATTENTION: works with hexahedrals only so far // calculate refined nodes positions (27), create nodes, return these nodes numbers

// linear <--> quadratic 
public:
	virtual void LinearToQuadratic(IVector& subset = IVector(0)); // make all linear elements quadratic - consistency !NOT! conserved unless all elements are linear in the beginning 
	virtual void QuadraticToLinear(IVector& subset = IVector(0),int removeunused = 1);	// make all quadratic elements linear - consistency !NOT! conserved unless all elements are quadratic in the beginning 
  int AddCenterNode(TArray<int>& points); // creates a node in the center of all nodes registered in array points, adds the index of this node to points
protected:
	virtual void QuadNodes(IVector& newpoints, FEElement* linelem); // returns nodelist for the quadratic element, adds intermediate nodes 
  virtual void LinearNodes(IVector& newpoints, FEElement* linelem); // returns nodelist for the linear element 
	int AddIntermediateNode(int i, int j); // creates intermediate node between two existing (global number) nodes

// automatic joints
public:
	virtual int JoinAtPlane(FEMesh_Set& master_elements, FEMesh_Set& slave_elements, Vector3D point, Vector3D normal); //join all nodes of slave to faces of master, returns number of AreaContact
	virtual int JoinBodiesAtPlane(int master_body, int slave_body, Vector3D point, Vector3D normal); //join all nodes of slave to faces of master, returns number of AreaContact
	virtual int JoinMaterialAtPlane(int master_material, int slave_material, Vector3D point, Vector3D normal); //join all nodes of slave to faces of master, returns number of AreaContact
//TODO: virtual int JoinSelectionsAtPlane(FEMesh_Set& master, FEMesh_Set& slave, Vector3D point, Vector3D normal);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ operations on sets and single entities
//elements
public: 
  virtual void SetElementDomain(FEMesh_Set& in_set_elements, int in_int_domain_number); // changes the domain of a set of elements, adds additional nodes, also updates faces if they are associated with elements
	virtual void SetElementDomain(int in_int_element_number, int in_int_domain_number) { SetElementDomain(FEMesh_Set(TSetElements,IntVec1(in_int_element_number)), in_int_domain_number); }
	virtual void SetElementColor(FEMesh_Set& in_set_elements, Vector3D in_V3D_new_color); // changes the elementcolor for a set of elements
	virtual void SetElementColor(int in_int_element_number, Vector3D in_V3D_new_color) { SetElementColor(FEMesh_Set(TSetElements,IntVec1(in_int_element_number)), in_V3D_new_color); }

	virtual void PrismTo3Hexes(int elemnr); // splits 1 prism into 3 hexehedrals - NOT IMPLEMENTED YET

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ dimensionality of node
public:
	virtual int Is3DNode(int nodenumber); // determine if the Node is used in 3D elements only
	virtual int Is2DNode(int nodenumber); // determine if the Node is used in 2D elements only
//+ true element type
public:
	virtual int ElementIsHex(int elemnr); // determines if the Finite Element element is really a hexahedral
	virtual int ElementIsPrism(int elemnr); // determines if the Finite Element element is a prism
	virtual int ElementIsPyram(int elemnr); // determines if the Finite Element element is a pyramid

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ convert sets functions - use FEMesh_Set only
protected:
	virtual void ComputeNodesOfElements(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_elements, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes) ); // computes a set of nodes for a given set of elements, subset of valid nodes can be included 
	virtual void ComputeElementsOfNodes(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements) ); // computes a set of ALL elements for a given set of nodes, subset of valid elements can be included 
	virtual void ComputeElementsFromNodes(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements) ); // computes a set of elements that can be built with a given set of nodes, subset of valid elements can be included 

	virtual void ComputeNodesOfFaces(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_faces, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes) ); // computes a set of nodes for a given set of faces, subset of valid nodes can be included
	virtual void ComputeFacesOfNodes(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_faces = FEMesh_Set(TSetFaces) ); // computes a set of ALL faces for a given set of nodes, subset of valid faces can be included
	virtual void ComputeFacesFromNodes(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_faces = FEMesh_Set(TSetFaces) ); // computes a set of faces that can be built a given set of nodes, subset of valid faces can be included

	virtual void ComputeFacesOfElements(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_elements, FEMesh_Set& subset_faces = FEMesh_Set(TSetFaces) ); // computes a set of faces for a given set of elements, subset of valid faces can be included
	virtual void ComputeElementsOfFaces(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_faces, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements) ); // computes a set of elements for a given set of faces, subset of valis elements can be included
	
	virtual void ComputeMaterialsOfElements(FEMesh_Set& rv_set_materials, FEMesh_Set& in_set_elements, FEMesh_Set& subset_materials = FEMesh_Set(TSetMaterials) ); // computes a set of all materials used in a given set of elements, subset of valid materials can be included
	virtual void ComputeElementsOfMaterials(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_materials, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements) ); // computes a set of all elements with materials in a given set, subset of valid elements can be included
	
	virtual void ComputeBodiesOfElements(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_elements, FEMesh_Set& subset_bodies = FEMesh_Set(TSetBodies) ); // computes a set of all bodies used in a given set of elements, subset of valid bodies can be included
	virtual void ComputeElementsOfBodies(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements) ); // computes a set of all elements with bodynumbers in a given set, subset of valid elements can be included
	
	virtual void ComputeBodiesOfNodes(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_bodies = FEMesh_Set(TSetBodies) ); // computes a set of all bodies used in a given set of nodes, subset of valid bodies can be included
	virtual void ComputeNodesOfBodies(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes) ); // computes a set of all nodes with bodynumbers in a given set of bodies

// not directly implemented, 2-stage via elements, no subset
	virtual void ComputeNodesOfMaterials(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_materials, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfMaterials(temp_elements, in_set_materials);
		ComputeNodesOfElements(rv_set_nodes, temp_elements, subset_nodes);
	}
	
	virtual void ComputeMaterialsOfNodes(FEMesh_Set& rv_set_materials, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_materials = FEMesh_Set(TSetMaterials) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfNodes(temp_elements, in_set_nodes);
		ComputeMaterialsOfElements(rv_set_materials, temp_elements, subset_materials);
	}

	virtual void ComputeFacesOfBodies(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_faces = FEMesh_Set(TSetFaces) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfBodies(temp_elements, in_set_bodies);
		ComputeFacesOfElements(rv_set_faces, temp_elements, subset_faces);
	}
	
	virtual void ComputeBodiesOfFaces(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_faces, FEMesh_Set& subset_bodies = FEMesh_Set(TSetBodies) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfFaces(temp_elements, in_set_faces);
		ComputeBodiesOfElements(rv_set_bodies, temp_elements, subset_bodies);
	}

	virtual void ComputeFacesOfMaterials(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_materials, FEMesh_Set& subset_faces = FEMesh_Set(TSetFaces) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfMaterials(temp_elements, in_set_materials);
		ComputeFacesOfElements(rv_set_faces, temp_elements, subset_faces);
	}

	virtual void ComputeMaterialsOfFaces(FEMesh_Set& rv_set_materials, FEMesh_Set& in_set_faces, FEMesh_Set& subset_materials = FEMesh_Set(TSetMaterials) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfFaces(temp_elements, in_set_faces);
		ComputeMaterialsOfElements(rv_set_materials, temp_elements, subset_materials);
	}

	virtual void ComputeBodiesOfMaterials(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_materials, FEMesh_Set& subset_bodies = FEMesh_Set(TSetBodies) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfMaterials(temp_elements, in_set_materials);
		ComputeBodiesOfElements(rv_set_bodies, temp_elements, subset_bodies);
	}

	virtual void ComputeMaterialsOfBodies(FEMesh_Set& rv_set_materials, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_materials = FEMesh_Set(TSetMaterials) )
	{
		FEMesh_Set temp_elements(TSetElements);
		ComputeElementsOfBodies(temp_elements, in_set_bodies);
		ComputeMaterialsOfElements(rv_set_materials, temp_elements, subset_materials);
	}

// convert a set of any type into a set of desired type, flag to force clean set otherwise old entries are kept
	virtual void ConvertSet(FEMesh_Set& the_set, TSetType typei, int flag_flushall = 0); 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ convert-access functions to be used in model files - allow several input/output entities - no calculation
public:
// Nodes:
	virtual void GetNodesOfElements(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_elements, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfElements(rv_set_nodes, in_set_elements, subset_nodes);
	}
	virtual void GetNodesOfElement(FEMesh_Set& rv_set_nodes, int in_int_elementnumber, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfElements(rv_set_nodes, FEMesh_Set(TSetElements, IntVec1(in_int_elementnumber)), subset_nodes);
	}
  
	virtual void GetNodesOfFaces(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_faces, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfFaces(rv_set_nodes, in_set_faces, subset_nodes);
	}
  virtual void GetNodesOfFaces(FEMesh_Set& rv_set_nodes, int in_int_facenumber, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfFaces(rv_set_nodes, FEMesh_Set(TSetFaces,IntVec1(in_int_facenumber)), subset_nodes);
	}
	virtual void GetNodesOfArea(FEMesh_Set& rv_set_nodes, int in_int_areanumber, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfFaces(rv_set_nodes, GetArea(in_int_areanumber), subset_nodes);
	}
	virtual void GetNodesOfSurface(FEMesh_Set& rv_set_nodes, FEMesh_Set subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfFaces(rv_set_nodes, Surface(), subset_nodes);
	}

	virtual void GetNodesOfMaterials(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_materials, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfMaterials(rv_set_nodes, in_set_materials, subset_nodes);
	}
	virtual void GetNodesOfMaterial(FEMesh_Set& rv_set_nodes, int in_int_materialnumber, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfMaterials(rv_set_nodes, FEMesh_Set(TSetMaterials, IntVec1(in_int_materialnumber)), subset_nodes);
	}

	virtual void GetNodesOfBodies(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfBodies(rv_set_nodes, in_set_bodies, subset_nodes);
	}
	virtual void GetNodesOfBody(FEMesh_Set& rv_set_nodes, int in_int_bodynumber, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes))
	{
		ComputeNodesOfBodies(rv_set_nodes, FEMesh_Set(TSetBodies, IntVec1(in_int_bodynumber)), subset_nodes);
	}

// Elements:
	virtual void GetElementsOfNodes(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements))
	{
		ComputeElementsOfNodes(rv_set_elements, in_set_nodes, subset_elements);
	}
// NO GETELEMENTSOFNODE HERE
	//virtual void GetElementsOfFaces();
	
	virtual void GetElementsOfMaterials(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_materials, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements))
	{
		ComputeElementsOfMaterials(rv_set_elements, in_set_materials, subset_elements);
	}
	virtual void GetElementsOfMaterial(FEMesh_Set& rv_set_elements, int in_int_materialnumber, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements))
	{
		ComputeElementsOfMaterials(rv_set_elements, FEMesh_Set(TSetMaterials,IntVec1(in_int_materialnumber)), subset_elements);
	}

	virtual void GetElementsOfBodies(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements))
	{
		ComputeElementsOfBodies(rv_set_elements, in_set_bodies, subset_elements);
	}
	virtual void GetElementsOfBody(FEMesh_Set& rv_set_elements, int in_int_bodynumber, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements))
	{
		ComputeElementsOfBodies(rv_set_elements, FEMesh_Set(TSetBodies, IntVec1(in_int_bodynumber)), subset_elements);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ geometry functions
// nodes
public:
	virtual int GetNodeAtPos(Vector3D pos, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL); //get (first) nodenumber for a given position, 0 if no node is found
	virtual int GetNodeAtPos(Vector2D pos, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL)	{ return GetNodeAtPos(pos.MakeV3D(), subset_nodes, tol); } //get (first) nodenumber for a given position

  virtual void GetNodesOnPlane(FEMesh_Set& rv_set_nodes, Vector3D nplane, double cplane, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL); // get nodes on a plane (as HNF-equation)
  virtual void GetNodesOnPlane(FEMesh_Set& rv_set_nodes, Vector2D nplane, double cplane, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL)	{	GetNodesOnPlane(rv_set_nodes, nplane.MakeV3D(), cplane, subset_nodes, tol);	} // get nodes within tolerance of a plane (as HNF-equation)
	virtual void GetNodesMinRespVect(FEMesh_Set& rv_set_nodes, Vector3D nplane, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL); //get a list of nodes that have minimal vaule of distance to plane (min d in HNF-equation)
	virtual void GetNodesMinRespVect(FEMesh_Set& rv_set_nodes, Vector2D nplane, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL) { GetNodesMinRespVect(rv_set_nodes, nplane.MakeV3D(), subset_nodes, tol); } //get a list of nodes that have minimal vaule of distance to plane (min d in HNF-equation)

	virtual int SortNodesByCoordinate(FEMesh_Set& set_to_be_sorted, int coordinate); // sort nodes by their position / coordinate (ascending)

	virtual void GetNodesInBox(FEMesh_Set& rv_set_nodes, Box3D& box, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL); // get all nodes in a box 
	virtual void GetNodesInBox(FEMesh_Set& rv_set_nodes, Vector3D corner1, Vector3D corner2, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL) { GetNodesInBox(rv_set_nodes, Box3D(corner1, corner2), subset_nodes, tol); } // get all nodes in a box 
	virtual void GetNodesInBox(FEMesh_Set& rv_set_nodes, Box2D& box, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL) { GetNodesInBox(rv_set_nodes, Box3D(box.PMin().MakeV3D(), box.PMax().MakeV3D()), subset_nodes, tol); }// get all nodes in a box 
	virtual void GetNodesInBox(FEMesh_Set& rv_set_nodes, Vector2D corner1, Vector2D corner2, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL) { GetNodesInBox(rv_set_nodes, Box3D(corner1.MakeV3D(), corner2.MakeV3D()), subset_nodes, tol); }// get all nodes in a box 

	virtual void GetNodesOnCircle(FEMesh_Set& rv_set_nodes_sorted, Vector3D center, Vector3D normal, double radius, int flag_DISC=0, double radius_hole=0., FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), int get_every_X_node=1, double tol = MESH_STD_TOL); // computes nodes that are on the circle/ on the disc		//$ DR 2012-08-20 radius_hole added
	virtual void SortNodesOnCircle(FEMesh_Set& rv_set_nodes, Vector3D center, Vector3D normal, double radius); // sorts the nodes by angle, first node defined 0°
	virtual void GetNodesOnCylinder(FEMesh_Set& rv_set_nodes_sorted, Vector3D center_bottom, Vector3D center_top, double radius, int flag_FULL=0, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), int get_every_X_node=1, double tol = MESH_STD_TOL); // computes nodes that are on the cylinder shell/ in the cylinder
	virtual void SortNodesOnCylinder(FEMesh_Set& rv_set_nodes, Vector3D center_bottom, Vector3D center_top, double radius); // sorts the nodes by angle and axial position

	virtual int GetNodesInCylinderShell(FEMesh_Set& rv_set_nodes, Vector3D base, Vector3D axis, double ri, double ro, double hb, double ht, FEMesh_Set& subset_nodes = FEMesh_Set(TSetNodes), double tol = MESH_STD_TOL); // identify all nodes in a cylindrical shell (node is in shell)
	virtual void ScaleMesh(double scaling_factor);	// coordinates of the nodes are scaled: new = scaling_factor*old

	// elements
public:
	virtual int GetElementsInBox(FEMesh_Set& rv_set_elements, Box3D& box, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements), double tol = MESH_STD_TOL); // identify all elements in a box
  virtual int GetElementsInCylinderShell(FEMesh_Set& rv_set_elements, Vector3D base, Vector3D axis, double ri, double ro, double hb, double ht, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements), double tol = MESH_STD_TOL); // identify all elements in a cylindrical shell (center point is in cylinder) 
	virtual int IsPointInCylinder(Vector3D point, Vector3D m1, Vector3D m2, double ri, double ro, double tol = MESH_STD_TOL); // checks if a point is in a cylinder shell

	virtual Box3D GetBox3DofElement(int elemnr); 	// returns a surrounding Box3D for a FElement
	virtual Box3D GetBox3DofElements(FEMesh_Set& in_set_elements); 	// returns a surrounding Box3D for a Set of FElement

// nodes and elements
public:
	virtual void GetElementsAndLocalNodesOfGlobalNode(FEMesh_Set& rv_set_elements, int in_int_globalnodenumber, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements)); // finds elements containing the given global node, also store the local node number in Nodes-Array of set
  virtual void GetElementsAndLocalCoordsOfGlobalNode(FEMesh_Set& rv_set_elements, TArray<Vector3D>& rv_array_localcoords, int in_int_globalnodenumber, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements)); // finds elements containing the given global node, also stores the local coordinates in Vector3D-Array
  virtual void GetElementsAndLocalCoordsOfGlobalPosition(FEMesh_Set& rv_set_elements, TArray<Vector3D>& rv_array_localcoords, Vector3D& in_V3D_globalposition, FEMesh_Set& subset_elements = FEMesh_Set(TSetElements)); // finds elements containing the given global node, also stores the local coordinates in Vector3D-Array

	virtual int2 GetFirstElementAndLocalNodeOfGlobalNode(int node);
	virtual int GetGlobalPositionsOfElementNodes(int elemnr, TArray<Vector3D>& points); // fills a list with the element nodes global positions
	virtual int GetElementAndLocalCoordOfGlobalPosition(Vector3D& global, Vector3D& local, IVector& subset = IVector()); // returns element number and local coordinate of a global position
	virtual Vector3D GetCenterOfElement(FEElement& el); // returns center of an element (global coords)
	virtual Vector3D GetCenterOfElement(int elem) { return GetCenterOfElement(GetElement(elem)); } // returns center of an element (global coords)
	virtual Vector3D GetLocalCoord(FEElement& el, Vector3D& globalcoord); // calculates local position in respect to a chosen element of a global vector 
	virtual Vector3D GetLocalCoord(int elem, Vector3D& globalcoord) { return GetLocalCoord(GetElement(elem),globalcoord); } // calculates local position in respect to a chosen element of a global vector 
  virtual Vector3D GetNodeLocalCoord(FEElement& el, int nodenumber_global); // calculates local position of a global node in the given element
	virtual Vector3D GetNodeLocalCoord(int elem, int nodenumber_global) { return GetNodeLocalCoord(GetElement(elem),nodenumber_global); } // calculates local position of a global node in the given element

// elements

// faces
public:
	virtual double GetFaceArea(FEMesh_Face& feface); // returns Area or Length of a Side of an Finite Element
	virtual double GetElementFaceArea(int elem,int side); // returns Area or Length of a Side of an Finite Element
	virtual double ComputeAngleAtNode(int nodenumber, int4& face); // calculate the interior angle at a node for a given face
	virtual double ComputeNodeArea(int nodenumber, int4& face); 	// Computes part of the faces ares that is associated with the node (4-gon: node, sidecenter, facecenter, sidecenter)
	virtual int GetFaceNodePositions(int4& face, Vector3D& p1, Vector3D& p2, Vector3D& p3,  Vector3D& p4); // returns number of valid entries in face and fills p1..p4 with nodepositions
	virtual Vector3D GetFaceCenterPoint(int4& face); // returns the center point of the face ( Center of Gravity )
	virtual Vector3D GetPrevNodePos(int thisnodenr, int4& face); // returns the position of the (cyclic) previous node
	virtual Vector3D GetNextNodePos(int thisnodenr, int4& face); // returns the position of the (cyclic) next node

// general - vector - distance
//
 	virtual double DistanceFromPoint(int nodenum,Vector3D refpoint) //distance of node i from point
	{ Vector3D dv; dv = refpoint-GetPoint3D(nodenum); return dv.Norm(); } 
	virtual double DistanceFromPoint(int nodenum,Vector2D refpoint) //distance of node i from point
	{ Vector2D dv; dv = refpoint-GetPoint2D(nodenum); return dv.Norm(); } 
	virtual double DistanceFromPlane(int nodenum,Vector3D nplane, double cplane) //distance of node i from plane
	{ return (nplane*GetPoint3D(nodenum)-cplane) / nplane.Norm();} 
	virtual double DistanceFromPlane(int nodenum,Vector2D nplane, double cplane) //distance of node i from plane
	{ return (nplane*GetPoint2D(nodenum)-cplane) / nplane.Norm();} 

	virtual double AxialDistance(Vector3D point, Vector3D fixed, Vector3D dir)
	{
		Vector3D pq = point - fixed;
		double factor = (pq * dir) / (dir*dir);
		Vector3D d = pq - (dir * factor);
		return d.Norm();
	}
	
	virtual Vector2D NormVector2D(Vector2D p1,Vector2D p2) // returns 2D normal vector of the vector between the 2 points
	{ return Vector2D( p2.Y()-p1.Y() , -(p2.X()-p1.X()) ); }
	virtual Vector2D NormVector2D(int2 globalnodenums) // returns 2D normal vector of the vector between the 2 points
	{ return NormVector2D(GetPoint2D(globalnodenums(2)),GetPoint2D(globalnodenums(2))); }
	virtual Vector2D NormVector2D(FEMesh_Face& face)
	{ return NormVector2D(int2(face.GetNode(1),face.GetNode(2))); }

	virtual Vector3D NormVector3D(Vector3D p1,Vector3D p2,Vector3D p3,Vector3D p4) // returns 3D normal vector the points or nullvector if they are not in one plane
	{
		if (IsPlanar(p1,p2,p3,p4)) return (p2-p1).Cross((p3-p1));
		else return Vector3D(0.0,0.0,0.0);
	}
	virtual Vector3D NormVector3D(int4 globalnodenums) // returns 3D normal vector the points or nullvector if they are not in one plane
	{	
		if (globalnodenums(4) == 0) globalnodenums(4) = globalnodenums(1); // in case of trig
		return NormVector3D(GetPoint3D(globalnodenums(1)),GetPoint3D(globalnodenums(2)),GetPoint3D(globalnodenums(3)),GetPoint3D(globalnodenums(4)));
	}
	virtual Vector3D NormVector3D(FEMesh_Face& face)
	{	return NormVector3D((int4)face.GetNodes()); }

	virtual int IsPlanar(Vector3D p1,Vector3D p2,Vector3D p3,Vector3D p4); // returns 1 if all points are in one plane 
	virtual int IsPlanar(int4 globalnodenums) // returns 1 if all points are in one plane 
	{
		if (globalnodenums(4) == 0) globalnodenums(4) = globalnodenums(1); // in case of trig
		return IsPlanar(GetPoint3D(globalnodenums(1)),GetPoint3D(globalnodenums(2)),GetPoint3D(globalnodenums(3)),GetPoint3D(globalnodenums(4)));
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ cleanup
	virtual void FindUnusedNodes(IVector& unusednodes = IVector()); // returns list of nodes new numbers (... a(n) = 0 : delete node)
	virtual int RemoveUnusedNodes(); // remove nodes that are not used in any element
	virtual void FindDoubleNodes(IVector& doublenodes = IVector(),double tol = 1E-10); // returns list of nodes new numbers (... a(n+1) <= a(n) : merge node)
	virtual int RemoveDoubleNodes(double tol = 1E-10); // merge nodes that occupy the same position (with tolerance)
	virtual int RemoveNodes(IVector& newnodenumebers); // actually removes the nodes, changes nodenumbers in elements accordingly

	virtual int RemoveUnusedFaces(); // removes faces with entries node == 0 from the list
	
	virtual void EraseNode(int i) { points.Erase(i); 	InitializeNodeSearchTree(); } // erases a node from nodelist list
	virtual int EraseManyNodes(IVector& newnodenumbers) { return points.RearrangeArray(newnodenumbers); InitializeNodeSearchTree(); } // erase many nodes from the dataarray (built array all new)
	virtual void EraseFace(int i); // erase face from the dataarray
	virtual void EraseElem(int i); // erase element from the dataarray
	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ auxiliary functions
//make functions for all FEElements

	FEElement* MakeTet(TArray<int> points, int material) //new Tet
	{
		FEElement* elem = new FETet(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeTetquad(TArray<int> points, int material) //new Tetquad
	{
		FEElement* elem = new FETetquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeHex(TArray<int> points, int material) //new Hex
	{
		FEElement* elem = new FEHex(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeHexquad(TArray<int> points, int material) //new Hexquad
	{
		FEElement* elem = new FEHexquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	// $EK 2013-03-11 added for prisms
	FEElement* MakePrism(TArray<int> points, int material) //new Prism
	{
		FEElement* elem = new FEPrism(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	// $EK 2013-03-05 added
	FEElement* MakePrismquad(TArray<int> points, int material) //new Prismquad
	{
		FEElement* elem = new FEPrismquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	//$ DR 2013-03-13 added
	FEElement* MakePyramid(TArray<int> points, int material) //new Pyramid
	{
		FEElement* elem = new FEPyramid(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	// $EK 2013-03-05 added
	FEElement* MakePyramidquad(TArray<int> points, int material) //new Pyramidquad
	{
		FEElement* elem = new FEPyramidquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeTrig(TArray<int> points, int material) //new Trig
	{
		FEElement* elem = new FETrig(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeTrigquad(TArray<int> points, int material) //new Trigquad
	{
		FEElement* elem = new FETrigquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeQuad(TArray<int> points, int material) //new Quad
	{
		FEElement* elem = new FEQuad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeQuadquad(TArray<int> points, int material) //new Quadquad
	{
		FEElement* elem = new FEQuadquad(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeLine(TArray<int> points, int material) //new Line
	{
		FEElement* elem = new FELine(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}
	FEElement* MakeLinequad(TArray<int> points, int material) //new Linequad
	{
		FEElement* elem = new FELine(points);
		elem->MaterialNum() = material;
		elem->Domain() = 1;
		return elem;
	}	

protected: // VARIABLES
// concerning MBS
	MBS * mbs;	                          // the MBS
	mutable FEMesh_Settings settings;   // settings used when converting the simpified mesh elements to MBS elements
	                                    // this member is declared as "mutable" to avoid changes in the client code
	TArray<int> mbsnodenumbers;					// stores mbs materialnumbers for the mesh-materials, created by AddMaterialsToMBS()
	TArray<int> mbsmaterialnumbers;			// stores mbs nodenumbers for the mesh-nodes, created by AddNodesToMBS()
	TArray<int> mbselementnumbers;      // stores mbs elementnumbers for the mesh-elements, created by AddElementsToMBS()

// mesh data
	TArray<FEMesh_Node> points;         // Nodes
	TArray<FEElement*> elements;				// Elements
	TArray<FEMesh_Face> faces;          // Faces (sides of single elements)
	
	TArrayDynamic<FEMesh_Set> areas;		// store a list of Areas, each is a set of faces
	FEMesh_Set surface;                 // determined to hold surface ( special set )

// materials
	TArray<Material*> material_data;    

// boundary conditions
  TArray<FEMesh_Load*> load_data;
	TArray<FEMesh_Constraint*> constraint_data;

// initial conditions
	TArray<Vector3D> init_vel;					//initial nodal velocities, initially zero length
 	TArray<Vector3D> init_disp;					//initial nodal displacements, initially zero length

// quick search and mapping
	Box3D nodebox;                      // box enveloping all nodes (must be initialized)
	SearchTree nodestree, add2mbstree;  // Searchtree for Mesh and additional tree for use in NodesToMBS routines ( no accessfunctions for that )
	TArray<IVector*> nodes_to_elements; // list elements for each node
	FEMesh_Set currentselection;        // buffer for all kind of selections selections 
	FEMesh_Set mesh_register;           // determined to hold lists of all { materials / domains / nodes } actually used by elements 
	                                    // and { faces } actually used by areas
	TArrayDynamic<FEMesh_Set> selections;      // list to save selections - searchable by name ( due to operator == of FEMesh_Set )	

// generator class
	FEMesh_Generator thegenerator;      // generator class for blocks and parts, generates the nodes elements in the parent mesh...


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE REMOVED  TO BE +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FUNCTION FORWARDS TO GENERATOR - THESE FUNCTIONS ARE CALLED IN MODEL***.CPP
// WHEN DIRECT CALLS ARE REPLACED BY INSTANCES OF BLOCKS THESE FORWARDS CAN BE FEMOVED
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		  IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																			MBSLoad& bodyload, int bodyloadactive,
		  																IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			int order=1)
	{ 
		mbs->UO(UO_LVL_warn) << "WARNING: GenerateHexahedralBlock(...) depreciated, please update your model file to 'Generator().GenerateHexahedralBlock(...)' (AD)\n";
		return Generator().GenerateHexahedralBlock(block, divisions_xyz, bodynr, matnr, color, ret_elementnumbers, ret_nodenumbers, ret_faces, constraintfaces_constainttypes,
	                                             constraintfaces_directions, bodyload, bodyloadactive, loadfaces_active, loadfaces_loadstrength, order);
	}
	virtual int GenerateDisc(double radius, double radius_hole, double height, int divisions_angle, int divisions_radial, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color,
							 IVector& blockelements, IVector& blocknodes)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: GenerateHexahedralBlock(...) depreciated, please update your model file to 'Generator().GenerateHexahedralBlock(...)' (AD)\n";
		return Generator().GenerateDisc(radius, radius_hole, height, divisions_angle, divisions_radial, divisions_height, pos, bodynr, matnr, color, blockelements, blocknodes);
	}


	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces)		
	{	
		mbs->UO(UO_LVL_warn) << "WARNING: GenerateHexahedralBlock(...) depreciated, please update your model file to 'Generator().GenerateHexahedralBlock(...)' (AD)\n";
		return Generator().GenerateHexahedralBlock(block, divisions_xyz, bodynr, matnr, color, blockelements, blocknodes, blockfaces);			
	}
	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: GenerateHexahedralBlock(...) depreciated, please update your model file to 'Generator().GenerateHexahedralBlock(...)' (AD)\n";
		return Generator().GenerateHexahedralBlock(block, divisions_xyz, bodynr, matnr, color);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FUNCTION FORWARDS TO FUNCTIONS THAT USE FEMESH_SET CLASS - THESE FUNCTIONS ARE CALLED IN MODEL***.CPP
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 // virtual void GetNodesOfElements(TArray<int>& elems, TArray<int>& nodes)
	//{ 
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOfElements(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set in_set_elements(TSetElements,elems), rv_set_nodes(TSetNodes);
	//	GetNodesOfElements(rv_set_nodes,in_set_elements);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesOfElements(int elemnr, TArray<int>& nodenumbers)
	//{
	//	GetNodesOfElements(IntVec1(elemnr), nodenumbers);
	//}

	//virtual void GetNodesOfFace(int facenumber, IVector& nodelist)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOfFace(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes);
	//	GetNodesOfFaces(rv_set_nodes,facenumber);
	//	nodelist.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesOfArea(int areanumber, IVector& nodelist)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOfArea(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes);
	//	GetNodesOfArea(rv_set_nodes,areanumber);
	//	nodelist.CopyFrom(rv_set_nodes.Nodes());
	//}

 // virtual void GetNodesOfBodies(TArray<int>& bodies, TArray<int>& nodelist)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOfBodies(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set in_set_bodies(TSetBodies,bodies), rv_set_nodes(TSetNodes);
	//	GetNodesOfBodies(rv_set_nodes,in_set_bodies);
	//	nodelist.CopyFrom(rv_set_nodes.Nodes());
	//}

 // virtual void GetElementswithMaterial(int material, TArray<int>& elements)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetElementswithMaterial(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_elements(TSetElements);
	//	GetElementsOfMaterial(rv_set_elements, material);
	//	elements.CopyFrom(rv_set_elements.Elements());
	//}

 // virtual void GetElementsOfBodies(TArray<int>& bodies, TArray<int>& elements)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetElementsOfBodies(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set in_set_bodies(TSetBodies, bodies),rv_set_elements(TSetElements);
	//	GetElementsOfBodies(rv_set_elements,in_set_bodies);
	//	elements.CopyFrom(rv_set_elements.Elements());
	//}

//+ geometry functions
	//virtual void GetNodesInBox(IVector& nodes, Box3D box, IVector& subset = IVector())
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesInBox(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesInBox(rv_set_nodes, box, subset_nodes);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

 // virtual void GetNodesInBox(IVector& nodes, Vector3D corner1, Vector3D corner2, IVector& subset = IVector())
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesInBox(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesInBox(rv_set_nodes, Box3D(corner1, corner2), subset_nodes);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesInBox(IVector& nodes, Box2D box, IVector& subset = IVector())
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesInBox(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	Box3D box3d(Vector3D(box.PMin().X(),box.PMin().Y(), MESH_STD_TOL), Vector3D(box.PMax().X(),box.PMax().Y(), MESH_STD_TOL)); 
	//	GetNodesInBox(rv_set_nodes, box3d, subset_nodes);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesInBox(IVector& nodes, Vector2D corner1, Vector2D corner2, IVector& subset = IVector())
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesInBox(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	Box3D box3d(Vector3D(corner1.X(),corner1.Y(), MESH_STD_TOL), Vector3D(corner2.X(),corner2.Y(), MESH_STD_TOL)); 
	//	GetNodesInBox(rv_set_nodes, box3d, subset_nodes);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}
 
	//virtual void GetNodesOnPlane(IVector& nodes, Vector3D nplane, double cplane, IVector& subset = IVector(), double tol = MESH_STD_TOL)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOnPlane(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesOnPlane(rv_set_nodes, nplane, cplane, subset_nodes, tol);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesOnPlane(IVector& nodes, Vector2D nplane, double cplane, IVector& subset = IVector(), double tol = MESH_STD_TOL)
	//{
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOnPlane(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesOnPlane(rv_set_nodes, nplane, cplane, subset_nodes, tol);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesMinRespVect(IVector& nodes, Vector3D nplane, IVector& subset = IVector(), double tol = MESH_STD_TOL)
 // {
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesMinRespVect(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesMinRespVect(rv_set_nodes, nplane, subset_nodes, tol);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

	//virtual void GetNodesMinRespVect(IVector& nodes, Vector2D nplane, IVector& subset = IVector(), double tol = MESH_STD_TOL)
 // {
	//	mbs->UO(UO_LVL_warn) << "WARNING: GetNodesMinRespVect(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
	//	FEMesh_Set rv_set_nodes(TSetNodes), subset_nodes(TSetNodes,subset);
	//	GetNodesMinRespVect(rv_set_nodes, nplane, subset_nodes, tol);
	//	nodes.CopyFrom(rv_set_nodes.Nodes());
	//}

public:
	//get a list of nodes on a a circle (center M, radius r, vector n) - flag_IN for entire area 	
	//if get_every_X_node is set to 1, all nodes are returned. otherwise every 2nd, 3rd, ... node
	virtual void GetNodesOnCircle(IVector& nodelist, Vector3D M, Vector3D n, double r, IVector& subset = IVector(), double tol = MESH_STD_TOL, int get_every_X_node=1, int flag_IN=0)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOnCircle(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
		FEMesh_Set rv_set_nodes(TSetNodes, nodelist);
		FEMesh_Set subset_nodes(TSetNodes, subset);
	  GetNodesOnCircle(rv_set_nodes, M, n, r, flag_IN, 0.,subset_nodes, get_every_X_node, tol);
		nodelist.CopyFrom(rv_set_nodes.Nodes());
	}
	//sort nodes in nodelist, it is assumed that they are on a circle with center M radius r and normalvector n
	//nodes are sorted according to the angle
	virtual void SortNodesOnCircle(IVector& nodelist, Vector3D M, Vector3D n, double r)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: SortNodesOnCircle(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
		FEMesh_Set rv_set_nodes(TSetNodes, nodelist);
		SortNodesOnCircle(rv_set_nodes, M, n, r);
		nodelist.CopyFrom(rv_set_nodes.Nodes());
	}

	//get a list of nodes on a a cylinder (center of circles M1 and M2, radius r)  - flag_IN for entire volume
	//if get_every_X_node is set to 1, all nodes are returned. otherwise every 2nd, 3rd, ... node
	virtual void GetNodesOnCylinder(IVector& nodelist, Vector3D M1, Vector3D M2, double r, IVector& subset = IVector(), double tol = MESH_STD_TOL, int get_every_X_node=1, int flag_IN=0)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: GetNodesOnCylinder(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
		FEMesh_Set rv_set_nodes(TSetNodes, nodelist);
		FEMesh_Set subset_nodes(TSetNodes, subset);
	  GetNodesOnCylinder(rv_set_nodes, M1, M2, r, flag_IN, subset_nodes, get_every_X_node, tol);
		nodelist.CopyFrom(rv_set_nodes.Nodes());
	}

	//sort nodes in nodelist, it is assumed that they are on a cylinder with centers of circles M1 and M2 and radius r 
	//nodes are sorted according to the angle and the axial position
	virtual void SortNodesOnCylinder(IVector& nodelist, Vector3D M1, Vector3D M2, double r)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: SortNodesOnCylinder(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
		FEMesh_Set rv_set_nodes(TSetNodes, nodelist);
		SortNodesOnCircle(rv_set_nodes, M1, M2, r);
		nodelist.CopyFrom(rv_set_nodes.Nodes());
	}

// elements
public:
  virtual int GetElementsInBox(IVector& elements, Box3D& box, IVector& subset = IVector(), double tol = 1E-10)
	{
		mbs->UO(UO_LVL_warn) << "WARNING: GetElementsInBox(...) call with IVector& depreciated, use call with FEMesh_Set instead! (AD)\n";
		FEMesh_Set rv_set_elements(TSetElements, elements);
		FEMesh_Set subset_elements(TSetElements, subset);
		GetElementsInBox(rv_set_elements, box, subset_elements, tol);
		elements.CopyFrom(rv_set_elements.Elements());
		return rv_set_elements.NElements();
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CLEANUP DONE UP TO THIS LINE - CONTINUE CLEANUP BELOW !!! 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ derived class FEMesh2D:   base class FEMesh 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//triangular finite element mesh, for storing and data processing
class FEMesh2D: public FEMesh
{
private:
	void Abstractor() {;};
public:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle
	//FEMesh2D():FEMesh() 
	//{
	//}
	FEMesh2D(const FEMesh2D& other)
	{
		CopyFrom(other);
	}
	virtual ~FEMesh2D()
	{
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle II
	FEMesh2D(MBS * mbsi):FEMesh(mbsi) 
	{
	}
	virtual void CopyFrom(const FEMesh2D& e)
	{
		FEMesh::CopyFrom(e);
		boundaryedges = e.boundaryedges; // specific FEMesh2D
	}
	virtual FEMesh2D* GetCopy() const
	{
		FEMesh2D* ed = new FEMesh2D(*this);
		return ed;
	}
	virtual void Reset()
	{
		FEMesh::Reset();
		boundaryedges.SetLen(0);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// remaining specific functions
	virtual int GetMeshDim() const { return 2; }
	virtual void GetElementPoints(int i, TArray<Vector2D>& elpoints) const
	{
		elpoints.SetLen(0);
		for (int j = 1; j <= elements(i)->NNodes(); j++)
		{
			elpoints.Add(points(elements(i)->GetNode(j)).GetCoords2D());
		}
	}

	virtual int NBoundaryEdges() const {return boundaryedges.Length();}
	virtual int2 GetBoundaryEdge(int i) const {return boundaryedges(i);}

	virtual void LoadNetgenMesh2D(const mystr& filename, int type, int order,int domain_number = 1); //type: 1==trig, 2==quad; order: 1=linear, 2=quadratic, domain_number: 0=domain number from file, 1,2,3,..= set whole mesh to specified domain
	virtual void ComputeBoundary();

private:
	TArray<int2> boundaryedges; //node1, node2
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ derived class FEMesh3D:   base class FEMesh 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//tetrahedral finite element mesh, for storing and data processing
class FEMesh3D: public FEMesh
{ 
private:
	void Abstractor() {;};
public:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle
	//FEMesh3D():FEMesh() 
	//{
	//}
	FEMesh3D(const FEMesh3D& other)
	{
		CopyFrom(other);
	}
	virtual ~FEMesh3D()
	{
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ lifecycle II
	FEMesh3D(MBS * mbsi):FEMesh(mbsi) 
	{
	}
	virtual void CopyFrom(const FEMesh3D& e)
	{
		FEMesh::CopyFrom(e);
		elementCG = e.elementCG;
		elementVol = e.elementVol;
	}
	virtual FEMesh3D* GetCopy() const
	{
		FEMesh3D* ed = new FEMesh3D(*this);
		return ed;
	}
	virtual void Reset()
	{
		FEMesh::Reset();
		elementCG.SetLen(0);
		elementVol.SetLen(0);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// remaining specific functions
	virtual int GetMeshDim() const { return 3; }

	//virtual const Vector3D& Point(int i) const {return points(i).GetCoords3D();} 
	virtual Vector3D& Point(int i) {return points(i).Coords3D();} //$ AD this function should be removed! needed for some models.... (models3.cpp,..
		
	virtual void GetElementPoints(int i, TArray<Vector3D>& elpoints) const //all points of element
	{
		elpoints.SetLen(0);
		for (int j = 1; j <= elements(i)->NNodes(); j++)
		{
			elpoints.Add(points(elements(i)->GetNode(j)).GetCoords3D());
		}
	}

	//type: 1==tet; order: 1=linear, 2=quadratic
	//netgenmode=1: netgen from export mesh (neutral), netgenmode=0: from ANSYS
	virtual void LoadNetgenMesh3D(const mystr& filename, int order, int material_num=0, int bodynum=1); 
//	virtual void LoadNetgenMesh3D(const mystr& filename, int type, int order,
//		 double rho0, double Em0, double nu0, int netgenmode); 

	// Load mesh from ANSYS node and element files
	//type: 1==tet; 2== hex; order: 1=linear, 2=quadratic
	// so far, only linear hexes are supported
	virtual void LoadAnsysNodesAndElements(const string& filename_nodes, const string& filename_elements, int type, int order, int bodynum);

	//load files generated with Tool "Gene.mac" from MN:
	virtual void LoadAnsysNodesAndElementsWithCGandMaterial(const mystr& filename_nodes, const mystr& filename_elements);

	virtual void AddPointsToReferenceFrame(int frameind, int bodyindex);

	// transform the whole mesh
	// x_new = translation + rotation * x_old
	// ATTENTION: only mesh point list is transformed, nodes already added to mbs are not changed
	virtual void Transform(const Vector3D& translation, const Matrix3D& rotation, const IVector& subset = IVector(0));

	// deform the whole Mesh, by nonlinear shearing
	// Coordinate deformed_coord of all points will be changed w.r.t. MathFunction final_shape
	virtual void Distort(int independent_coord, int deformed_coord, TArray<MathFunction*> final_shape, TArray<double> ref_values, int deformed_coord2 = 0, const IVector& subset_nodes = IVector(0));

	virtual Vector3D GetElementGlobalCG(int elem) const {return elementCG(elem);}

	virtual double GetElementVolume(int elem) const {return elementVol(elem);}

	virtual void DrawElements(const TArray<Vector3D>& node_disp);

private:
//from special load routine 
	TArray<Vector3D> elementCG;  //center of gravity of element, with respect to global coordinates ==> used for CMS/FFRF method
	TArray<double> elementVol;	 //volume of each element ==> used for CMS/FFRF method

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ THE SCRAP YARD - these functions sould be eliminated a.s.a.p.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class FEPRISM: does no longer work with current base class / did never work anyway
//////class FEPrismquad:public FEElement //according to defs of ANSYS
//////{
//////public:
//////	enum { NNODES = 15 };
//////	enum { NFACES = 5 };
//////	enum { NFACES3 = 2 };
//////	enum { NFACES4 = 3 };
//////
//////	FEPrismquad():FEElement() {};
//////	FEPrismquad(TArray<int>& points):FEElement() 
//////	{
//////		for (int i=0; i<NNodes(); i++)
//////			node[i] = points(i+1);
//////	};
//////
//////	virtual FEPrismquad* GetCopy() const
//////	{
//////		FEPrismquad* ed = new FEPrismquad(*this);
//////		return ed;
//////	}
//////	FEPrismquad(const FEElement& e) {CopyFrom(e);}
//////	virtual void CopyFrom(const FEPrismquad& e)
//////	{
//////		FEElement::CopyFrom(e);
//////		for (int i=0; i<NNodes(); i++)
//////			node[i] = e.node[i];
//////	}
//////
//////	virtual int Type() const {return 4;} //1==trig, 2==quad, 3=tet, 4=hex, 5=prism
//////	virtual int Order() const {return 2;} //1==linear, 2==quadratic
//////	virtual int NNodes() const {return NNODES;}
//////	virtual int GetNode(int i) const {return node[(i-1)%NNODES];}
//////	virtual void SetNode(int i, int pnum) {node[(i-1)%NNODES] = pnum;}
//////	virtual int NSides3() const {return NFACES3;}
//////	virtual int NSides4() const {return NFACES4;}
//////	virtual int NSides() const {return NFACES;}
//////	virtual int3 GetSide3(int i) const //maybe extend for all 16 triangulated surfaces?
//////	{
//////		switch(i)
//////		{
//////		case 1: return int3(1,3,2); break;
//////		case 2: return int3(1,2,4); break;
//////		default: assert(0); return int3(0,0,0);
//////		}
//////	}
//////	virtual int4 GetSide4(int i) const 
//////	{
//////		switch(i)
//////		{
//////		case 1: return int4(1,5,7,3); break;
//////		case 2: return int4(2,4,8,6); break;
//////		case 3: return int4(1,2,6,5); break;
//////		default: assert(0); return int4(0,0,0,0);
//////		}
//////	}
//////	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}
//////	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
//////private:
//////	int node[NNODES];
//////};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions used before class meshedpart and segmented block were written - 
// !!!FIND AND REPLACE IN MODEL FILES !!!
//void CutBlocks(TArray<SegmentedHexBlock*>& blocks, TArray<double> xcut = 0, TArray<double> ycut = 0, TArray<double> zcut = 0); // cuts a set of blocks - compute cutting planes, check priority for occupation, fill data (cutting planes, populaiton array) in HexBlock-class
//void FindAllCuttingPlanes(TArray<Box3D>& boxes, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut); // returns sorted arrays of cutting planes with all box borders and all prefilled entries within limits of overall box
//void ComputeOccupationArray(TArray<Box3D>& boxes, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut, TArray<int>& occupation); // computes field which block occupies the volume (the latter block always gets priority)
//void SetBlockData(TArray<SegmentedHexBlock*>& blocks, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut, TArray<int>& occupation); // writes computed sectioning to the blocks

#endif // FEMESH__H