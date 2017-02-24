//#**************************************************************
//#
//# filename:             models_include.cpp
//#
//# author:               Johannes Gerstmayr
//#
//# generated:						Sept 2010
//# description:          Compile Models files from ModelsLib to "models_include.o"
//												(see "..\ModelsLib\all_models.h")
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
//#**************************************************************

#include "mbs_interface.h"
#include "myfile.h"
#include "femathHelperFunctions.h"

#include "geomelements.h"     
#include "csg_geometry.h"  
#include "MBSLoad.h"

#include "solversettings_auto.h"
#include "options_class_auto.h" //JG 2013-01-10
#include "elementdataaccess.h" //JG 2013-01-10


#include "element.h"
#include "node.h"
#include "../femesh/femesh.h"
#include "../femesh/femesh_aux.h"
#include "sensors.h"
#include "sensorsSpecific.h"
#include "../ServiceObjectsLib/sensorProcessorsSpecific.h"
#include "material.h"
#include "constraint.h"
#include "control.h"

#include "KinematicPairs.h"
#include "RigidBodyJoints.h"
#include "SpecialConstraints.h"

#include "SDA_constraints.h"
#include "contact2D.h"
#include "contact3D.h"
#include "rigid2D.h"
#include "rigid3D.h"
#include "control.h"

#include "Sim2Hotint.h"
//#include "../MBSKernelLib/optimization.h"
//#include "sensitivity.h"
//#include "performcomputation.h"

#include "Rigid3DKardan.h"
#include "referenceframe2D.h"
#include "referenceframe3D.h"
#include "CMSElement2D.h"
#include "BaseCMSElement.h"
#include "CMSElement.h"
#include "FiniteElement2D.h"
#include "FiniteElement3D.h"
#include "GCMSElement.h"
#include "FiniteElement3DFFRF.h"

//#include "ANCFBeam2D.h"
//#include "ANCFBeam3D.h"
#include "ANCFCable2D.h"
#include "ANCFCable3D.h"
#include "ANCFBeamShear2D.h"
#include "ANCFBeamShear3D.h"
#include "ANCFBeam3DTorsion.h"
#include "ANCFBeamShearFE2D.h"
//#include "ANCFBeamBE2D.h"

#include "ANCFplate3D.h"

#include "Beam2DFFRF.h"
#include "Beam2DaFFRF.h"
#include "Beam3D.h"

#include "Plate2D.h"
#include "FE3DHexTet.h"
#include "FE2DTriQuad.h"
#include "ANCFThinPlate3D.h"
#include "Beam2DFFRF.h"
//#include "ANCFPipe2D.h" //extended elements
#include "ANCFAxMovBeam2D.h"
#include "Truss3D.h"

//#include "element.h"
//#include "ANCFBeam3DTorsion.h"
#include "AverageConstraint.h"
#include "mass1D.h"
#include "GCMSRotorElement.h"




#include "..\ModelsLib\all_models.h"
