//#**************************************************************
//# filename:            BaseCMSElement.h
//#
//# author:              Gerstmayr, Vetyukov 
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
 

#ifndef BaseCMSElement3D__H
#define BaseCMSElement3D__H

typedef enum { TCMSStaticModeTranslation = 1, TCMSStaticModeRotation = 2  } TCMSStaticModeType;

class SolverParameters
{
public:
	SolverParameters()
	{
		DefaultInitialize();
	}

	SolverParameters(int solvertypei, double tolerancei, int maxiterationsi)
	{
		solvertype = solvertypei; tolerance = tolerancei; maxiterations = maxiterationsi; lambda_precond = 0; use_precond = 0;
	}

	void DefaultInitialize()
	{
		solvertype = 2;      // use LOBPCG solver
		tolerance = 1e-8;
		maxiterations = 10000;
		lambda_precond = 0;
		use_precond = 0;
	}

	void CopyFrom(const SolverParameters& sp)
	{
			solvertype = sp.solvertype;
			tolerance = sp.tolerance;
			maxiterations = sp.maxiterations;
			lambda_precond = sp.lambda_precond;
			use_precond = sp.use_precond;
	}

	int solvertype;    // 0 = direct HOTINT solver, 1 = Matlab solver (iterative)
	double tolerance;  // error tolerance for Matlab
	int maxiterations; // maximum number of iterations for Matlab
	double lambda_precond; // lobpcg: use preconditioner K + lambda M
	int use_precond; // use preconditioner flag
};

//this element contains the whole structure of the body (finite element numbers, nodes) in this class.
//all nodes have local indices in the CMSElement
template<class RIGID>
class BaseCMSElement: public ReferenceFrame3D<RIGID>
{
public:
	enum {CMSmaxDOF = 100}; //100 original
	enum {CMSmin_sparse_size = 72}; //could be determined from maximum element size ...


	// Constructor, TCMSflag is set
	BaseCMSElement(MBS* mbsi):ReferenceFrame3D<RIGID>(mbsi), 
		boundarynode(), boundarydofelem(), boundarydof(), nbmodes(0), nimodes(0), staticmode_nodes(0), staticmode_type(0), staticmode_loccoords(0),
		massnodes(0), massvalues(0), n_zeromodes(-1), Dr(0), read_eigenmodes_from_file(0)
	{
		type = (TMBSElement)(type|TCMSflag);
		solverparameters.DefaultInitialize();
		EVfile = mystr("");
		animation_mode = 0;
	};

	// Copy constructor
	BaseCMSElement(const BaseCMSElement& e):ReferenceFrame3D<RIGID>(e.mbs),
		boundarynode(), boundarydofelem(), boundarydof(), read_eigenmodes_from_file(0)  {CopyFrom(e);};

	// destructor
	~BaseCMSElement()
	{
		for (int i=1; i<=staticmode_nodes.Length(); i++)
			delete staticmode_nodes(i);
		for (int i=1; i<=staticmode_loccoords.Length(); i++)
			delete staticmode_loccoords(i);
	}

	// number of dofs for the unreduced system
	virtual int SOSFull() const 
	{
		int sosfull = 0;
		//nodal dof:
		for (int i=1; i <= nodes.Length(); i++) sosfull+=nodes(i)->SOS();

		//element dof:
		for (int i=1; i <= NFFRFElements(); i++) sosfull+=GetFFRFElement(i).SOSowned();

		return sosfull;
	}

	virtual const double& GetXactFull(int i) const;
	virtual double& GetXactFull(int i); //should not be used ...
	virtual const double& GetDrawValueFull(int i) const;

	virtual void SetNZeromodes(int nz) {n_zeromodes = nz;}

	// Add elements and nodes from a FEMeshInterface to the CMS element
	void AddFEMesh(const FEMeshInterface& femesh);

	virtual Element* GetCopy()
	{
		Element* ec = new BaseCMSElement(*this);
		return ec;
	}
	
	virtual void CopyFrom(const Element& e)
	{
		ReferenceFrame3D<RIGID>::CopyFrom(e);
		const BaseCMSElement& ce = (const BaseCMSElement&)e;

		nbmodes = ce.nbmodes;
		nimodes = ce.nimodes;

		n_zeromodes = ce.n_zeromodes;

		boundarydofelem.SetLen(0);
		for (int i=1; i<=ce.boundarydofelem.Length(); i++)
		{
			boundarydofelem.Add(ce.boundarydofelem(i));
		}	

		boundarydof.SetLen(0);
		for (int i=1; i<=ce.boundarydof.Length(); i++)
		{
			boundarydof.Add(ce.boundarydof(i));
		}	

		boundarynode.SetLen(0);
		for (int i=1; i<=ce.boundarynode.Length(); i++)
		{
			boundarynode.Add(ce.boundarynode(i));
		}	

		M = ce.M;
		K = ce.K;
		Msparse = ce.Msparse;
		Ksparse = ce.Ksparse;
		Mr = ce.Mr;
		Kr = ce.Kr;
		Hr = ce.Hr;
		intrhodudqr = ce.intrhodudqr;

		Dr = ce.Dr;

		Phi_CB = ce.Phi_CB;
		volume = ce.volume;


		staticmode_nodes.SetLen(ce.staticmode_nodes.Length());
		for (int i=1; i<=staticmode_nodes.Length(); i++)
			staticmode_nodes(i) = new TArray<int>(*ce.staticmode_nodes(i));
		staticmode_loccoords.SetLen(ce.staticmode_loccoords.Length());
		for (int i=1; i<=staticmode_loccoords.Length(); i++)
			staticmode_loccoords(i) = new TArray<Vector3D>(*ce.staticmode_loccoords(i));
		staticmode_type = ce.staticmode_type;

		massnodes = ce.massnodes;
		massvalues = ce.massvalues;

		solverparameters.CopyFrom(ce.solverparameters);
		EVfile = ce.EVfile;
		read_eigenmodes_from_file = ce.read_eigenmodes_from_file;
		animation_mode = ce.animation_mode;
	}

	// initialization of  FFRF matrices is done here
	virtual void Initialize();
	virtual int PerformNodeCheck() const {return 0;}
	virtual const char* GetElementSpec() const {return "BaseCMSElement";}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//must be called before Assemble!!!!!
	////////link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	//////virtual void DoModalAnalysis(const TArray<int2>& fixednodes);
	// link elements and compute ltg map for element -> xgfull
	virtual void ComputeElementLTG();
	// compute entries of list boundarydof (0 for internal, 1 for boundarydof, 2 for staticmodedof, 3 for fixeddof)
	virtual void ComputeBoundaryDofIndexList(const TArray<int2>& fixednodes, int sosfull);
	// resort own ltg map according to (B) and (I) dofs
	virtual void ResortLTG(int sosfull);
	// compute number of zeromodes from fixednodes
	virtual void ComputeNZeromodes();
	// compute/read eigenmodes, build classic Craig-Bampton matrix Phi_CB which contains modes columnwise
	// nb_modesclassic .. number of classic boundary modes attached to one dof
	// nb_modesstatic ... number of additional static modes
	// sosfull .......... number of unreduced dofs
	// nb_dofstotal ..... number of unreduced boundary dofs
	// nb_dofsstatic .... number of unreduced dofs which belong to additional static modes
	virtual void BuildPhiCB_CMS(Matrix& Phi_CB, int nb_modesclassic, int nb_modesstatic, int sosfull, int nb_dofstotal, int nb_dofsstatic);
	// compute full mass and stiffness matrices
	virtual void ComputeMassAndStiffnessMatrix(SparseMatrix& Msparse, SparseMatrix& Ksparse);
	// Eigenmode Computation in Matlab
	virtual void Compute_Phi_X_Matlab(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, Matrix& X, Vector& eigval, Matrix& eigmodes);
	// Eigenmode Computation with own direct solver in HOTINT
	virtual void Compute_Phi_X(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, Matrix& X, Vector& eigval, Matrix& eigmodes);
	// Eigenmode Computation with own iterative sparse solver in HOTINT
	virtual void Compute_Phi_X_Sparse(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, Matrix& X, Vector& eigval, Matrix& eigmodes);
	// Normalization of the computed Eigenvectors, normalization mode = 0 (max(abs(v))=1) or 1 (v'*v = 1)
	virtual void NormalizeEigenvectors(int normalization_mode, int nf, int nr, Matrix& eigmodes);
	// method overwritten for (G)CMSElement
	virtual void DoModalAnalysis(const TArray<int2>& fixednodes)
	{
		UO() << "Error in BaseCMSElement::DoModalAnalysis called for base class!\n";
	}

 //precompute mass and store
	virtual void ComputeMass();
		

	virtual const Element& GetFFRFElement(int i) const {return GetMBS()->GetAuxElement(FFRFelements(i));}
	virtual Element& GetFFRFElement(int i) {return GetMBS()->GetAuxElement(FFRFelements(i));}

	// number of boundary modes = classical boundary modes + additional static modes
	virtual int NBModes() const {return nbmodes;} //boundary modes
	// number of internal modes due to modal analysis
	virtual int NIModes() const {return nimodes;} //internal modes
	virtual void SetNIModes(int i) {nimodes = i;} //internal modes
	// total number of modes
	virtual int NModes() const {return NBModes()+NIModes();}

	// get number of zero modes, which are not used in CMS method
	virtual int NZeroModes() const {return n_zeromodes;}


	virtual int IsCMS() const {return 1;}
	virtual int IsGCMS() const {return 0;}

	// Add boundary dof number locelemdof which is local to Element elem
	virtual void AddBoundaryElemDOF(int elem, int locelemdof) {boundarydofelem.Add(int2(elem,locelemdof));}
	// Add all dofs of Node node as boundary dofs
	virtual void AddBoundaryNode(int node) 
	{
		boundarynode.Add(int2(node,1));
		boundarynode.Add(int2(node,2));
		boundarynode.Add(int2(node,3));
	}
	// Add all dofs of Node node as boundary dofs
	virtual void AddBoundaryNodeDOF(int node, int nodedof) 
	{
		boundarynode.Add(int2(node,nodedof));
	}


	// additional static mode,
	// rigid body translation or rotation prescribed for nodes from nodeset (type = TCMSStaticModeTranslation or TCMSStaticModeRotation)
	// case Translation: direction is the direction of translation
	// case Rotation:    direction is the axis of rotation, center is a point on the axis of rotation
	virtual void AddStaticMode(TArray<int>& nodeset, TCMSStaticModeType type, Vector3D direction, Vector3D center = Vector3D(0.,0.,0.))
	{
		staticmode_type.Add(type);
		staticmode_nodes.Add(new TArray<int>(nodeset));
		TArray<Vector3D> *loccoords = new TArray<Vector3D>;
		loccoords->SetLen(0);
		loccoords->Add(direction);
		if (center.Norm())
			loccoords->Add(center);
		staticmode_loccoords.Add(loccoords);
	}

	//virtual Vector3D GetCenter(const TArray<int> nodenums) const;

	// Set Eigenvalue solver
	virtual void SetUseMatlab() {solverparameters.solvertype = 1;}
	virtual void SetUseDirect() {solverparameters.solvertype = 0;}
	virtual void SetUseLOBPCG() {solverparameters.solvertype = 2;}
	virtual void SetEVSolver(int i) {solverparameters.solvertype = i;}
	virtual void SetSolverParameters (int solvertype, double tolerance, int maxiterations, int use_precond=0, double lambda_precond = 0)
	{
		solverparameters.solvertype = solvertype;
		solverparameters.tolerance = tolerance;
		solverparameters.maxiterations = maxiterations;
		solverparameters.lambda_precond = lambda_precond;
		solverparameters.use_precond = use_precond;
	}
	virtual void SetEigenmodesFromFile(mystr modefile)
	{
		EVfile = modefile;
	}
	virtual int UseEigenmodeFile() { if (EVfile == mystr("")) {return 0;} else {return 1;} }

	// Add additional mass m to node nodenum
	virtual void AddNodeMass(int nodenum, double m)
	{
		massnodes.Add(nodenum);
		massvalues.Add(m);
	}

	virtual void SetDiagonalInternalDamping(const Vector& DrI) {Dr = DrI;}
	virtual const Vector& GetDiagonalInternalDamping() {return Dr;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	// AD: Reads Eigenmode Data from File this->EVFile, Data can bewritten to this->Phi_CB and/or DataManager
	// use ReadEigenModesFile(1) in Initialize() (optional)
	// use ReadEigenModesFile(0,1) to write Eigenvectors to Datamanager ( in ModelsFile )
	//void ReadEigenModesFile(int flag_toPhiCB = 1, int flag_toDataManager = 0);
	// AP: This is not good, since * maybe we want to read data into some other matrix my_local_Phi_CB
	//                                 * the Phi_CB-matrix is not loaded into the Data manager, 
	//                                   but essentially unit vectors with 1 in the j-th component correspond to mode j
	// Two routines:
	// ReadEigenModesFile reads eigenmodes from file columnwise into eigenmodemat
	virtual void ReadEigenModesFile(Matrix& eigenmodemat, int& successful);
	// Load solution vectors representing different eigenmodes into Data Manager
	virtual void EigenModesToDataManager();
	// switch to animation-mode (no computation of matrices, eigenmodes, etc)
	virtual void SetAnimationMode(int flag){animation_mode = flag;}

		virtual Box3D GetElementBox() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBox());
		}
		return box;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBoxD());
		}
		return box;
	}


protected:

	//CMS:
	int nbmodes;
	int nimodes;

	int n_zeromodes;  // number of zero Eigenmodes, which are not used in CMS method (e.g. when free-free modes are used, n_zeromodes = 6

	Matrix Phi_CB; //transformation matrix reduced to full dof
	TArray<int2> boundarydofelem; //element number and local dof which are fixed
	TArray<int> boundarydof; //mark dof with 1 which are fixed
	TArray<int2> boundarynode;     //node number and local dof which are fixed

	Matrix M, K, Mr, Kr, Hr; //full and reduced system matrices
	Matrix intrhodudqr;      // reduced system matrix int rho du/dq for gravity load
	Vector Dr; //diagonal modal damping terms
	SparseMatrix Msparse, Ksparse;
	double volume; //for intdudq

	// staticmode_nodes contains all nodes for which a constraint of type "staticmodes_type" is applied
	// these nodes are then also not included in the interior (I) nodes for which modes are computed
	TArray<TCMSStaticModeType> staticmode_type;
	TArray<TArray<int>*> staticmode_nodes;
	TArray<TArray<Vector3D>*> staticmode_loccoords;

	// additional node mass:
	TArray<int> massnodes;      // node numbers to which an additional mass is added
	TArray<double> massvalues;     // value of mass which is added to respective node from massnodes


	// Eigenvalue solver options
	SolverParameters solverparameters;
	// File from which Eigenmodes are read if there/to which Eigenmodes are written if not there
	mystr EVfile;
	// 1 if Eigenmodes are available in eigenmode file, 0 else
	int read_eigenmodes_from_file;
	// 1 if animation_mode (no computation of FFRF-Matrices, Eigenmodes,...), 0 else
	int animation_mode;
};

#endif // BaseCMSElement3D__H