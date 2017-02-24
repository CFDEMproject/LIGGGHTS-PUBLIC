/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: workspace.h                                             *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/


#ifndef WORKSPACE_H
#define WORKSPACE_H

#include "matrices.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>   
#include <iomanip>  
#include <vector>


class System;
class Solver;

//The Struct holding all the data, as well as the solver and integrator ids
struct SysData{
	System *    system;         // a pointer to the actual multibody system
	int         solver;         // id of the solver to be used
	int         integrator;     // id of the integrator to be used (NOT used in the code)
};

class Workspace {

    //Private data/structs
	SysData * systemContainer;   // a pointer to the structs that host the multibody systems data
	int currentIndex;            // the current id of the (last) multibody system = number of multibody systems that are computed
	int maxAlloc;                // just for saving the number of currently allocated multibody systems
	bool mydebug;
	
public:
     Workspace();
     ~Workspace();
     
     //Public data/structs
     double Thalf;              // half time step size
     double Tfull;              // full time step size
     double ConFac;             // ??
     double KE_val;             // kinetic energy of the system
     int    FirstTime;          // indicator used for first-time setting of kinetic energy
     
     //Pubic member functions
     bool LoadFile(char* filename);
     
     bool SaveFile(char* filename, int index = -1);

     System* GetSystem(int index = -1);
     
     void AddSolver(Solver* s, int index = -1);        
    

     void LobattoOne(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space);	
     
     void LobattoTwo(double **&vcm,double **&omega,double **&torque, double **&fcm);	 
     
 
     bool MakeSystem(int& nbody, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space, int &njoint, int **&jointbody, double **&xjoint, int& nfree, int*freelist, double dthalf, double dtv, double tempcon, double KE, int flag);
     																																							
	  
	  bool SaveSystem(int& nbody, double *&masstotal, double **&inertia, double **&xcm, double **&xjoint, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space, double **&acm, double **&alpha, double **&torque, double **&fcm, int **&jointbody, int &njoint);
	  																		
	 bool MakeDegenerateSystem(int& nfree, int*freelist, double *&masstotal, double **&inertia, double **&xcm, double **&vcm, double **&omega, double **&ex_space, double **&ey_space, double **&ez_space);
     int getNumberOfSystems();
     
     void SetLammpsValues(double dtv, double dthalf, double tempcon);
     void SetKE(int temp, double SysKE);
	  
	  void RKStep(double **&xcm, double **&vcm,double **&omega,double **&torque, double **&fcm, double **&ex_space, double **&ey_space, double **&ez_space);	
	  
	  void WriteFile(char* filename);

private:
	void allocateNewSystem(); //helper function to handle vector resizing and such for the array of system pointers
};

#endif
