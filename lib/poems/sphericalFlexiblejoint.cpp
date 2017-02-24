/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: sphericalFlexiblejoint.cpp                                      *
 *      AUTHORS: Mingqiu Wu, Stefan Radl                                           * 
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
 

#include "sphericalFlexiblejoint.h"
#include "point.h"
#include "matrixfun.h"
#include "body.h"
#include "fastmatrixops.h"
#include "norm.h"
#include "eulerparameters.h"
#include "matrices.h"
#include <iomanip>
#include <math.h> 
#include <cmath>
#include <unistd.h>     //for using the function sleep    
#include <stdio.h>  


SphericalFlexiblejoint::SphericalFlexiblejoint(){
} 
SphericalFlexiblejoint::~SphericalFlexiblejoint(){
}

void SphericalFlexiblejoint::ForwardKinematics()
{
    Vect3 result1,result2,result3,result4;
    //call base class function
    SphericalJoint::ForwardKinematics();
    //call additional calculation
    
    BendingAndTwistingTorque(body1->n_C_k,body2->n_C_k,result1,result3);
    //Assign torque to current body and neighbor body, opposite sign
    result2 = -result1;
    FastAssign(result1,body1->btorque);
    FastAssign(result2,body2->btorque);
    result4 = -result3;
    FastAssign(result3,body1->ttorque);
    FastAssign(result4,body2->ttorque);
    //cout<<" what is Twisttorque==  "<< result1(1) << endl; 
    // printf ("outside function: %g \n", result1(2));
}

void SphericalFlexiblejoint::BendingAndTwistingTorque(Mat3x3 temp1, Mat3x3 temp2, Vect3 &btorque,Vect3 &ttorque)
{
   Vect3 resultA,resultB,result1,bangle;
   Vect3 resultC,resultD,result2,tangle;
   resultA.Zeros(),resultB.Zeros(),result1.Zeros(),bangle.Zeros();
   resultC.Zeros(),resultD.Zeros(),result2.Zeros(),tangle.Zeros();
   double norm_e=0,norm_a=0,norm_b=0,BangTemp=0,abdot=0;
   double norm_y=0,norm_c=0,norm_d=0,TangTemp=0,cddot=0;
   double Kb,Bendinertia,Twistinertia, Kt;
   const double L  = 1.0;
   const double Ey = 1e3;//5e3;
   const double G  = 1e4;
   const double PI = 3.141592653589793238463;
   const double radiusA  = 0.05;
   const double radiusB  = 0.05;
   
       //calculate orientation for fibres (axial direction for bending and vertical direction for twisting)
       FastMult(body1->n_C_k,(body1->GetPoint(2))->position,resultA); 
       FastNegMult(body2->n_C_k,(body2->GetPoint(1))->position,resultB);
       //resultA(1)=0.933;resultA(2)=0;resultA(3)=0.25;
       for(int m = 0; m < 3; m++)
       {  
          //resultA(m+1)  = temp1(m+1,1);	//orientation body1 y-direction
          //resultB(m+1)  = temp2(m+1,1);	//orientation body1 y-direction
          resultC(m+1)  = temp1(m+1,2);	//orientation body1 y-direction
          resultD(m+1)  = temp2(m+1,2);	//orientation body2 y-direction
       }	
       //cout<<" what is resultB==  "<< resultB(1) << "..."<<resultB(2)<< "..."<<resultB(3)<<endl; 
       //cout<<" what is resultA==  "<< resultA(1) << "..."<<resultA(2)<< "..."<<resultA(3)<<endl;	
       // cross product of the orientation vector between ibody and ibody+1
       FastCross(resultA,resultB,result1);        
       norm_e = Magnitude(result1);
       FastCross(resultC,resultD,result2);
       norm_y = Magnitude(result2); 
       //normolize the cross product.
       if(norm_e>1e-30){
          for(int n = 0; n < 3; n++) result1(n+1)  /=  norm_e; 
       }
       else{
          for(int n = 0; n < 3; n++) result1(n+1)   =  0;
       }
       if(norm_y>1e-50){
          for(int n = 0; n < 3; n++) result2(n+1)  /=  norm_y; 
       }
       else{
          for(int n = 0; n < 3; n++) result2(n+1)   =  0;
       }
       //dot product of the orientation vector between ibody and ibody+1
       for (int i=1;i<=3;i++) abdot   +=   resultA(i)*resultB(i);
       for (int i=1;i<=3;i++) cddot   +=   resultC(i)*resultD(i);
       norm_a  =  Magnitude(resultA);
       norm_b  =  Magnitude(resultB);
       norm_c  =  Magnitude(resultC);
       norm_d  =  Magnitude(resultD);
       if(norm_a*norm_b>1e-30)
       {
              BangTemp      =  acos(abdot /(norm_a*norm_b));
              //cout<<" what is BangTemp==  "<< BangTemp << endl; 
       }
       else
       {
              BangTemp      =  0;
       }
       
       if(norm_c*norm_d>1e-50)
       {
              TangTemp      =  acos(cddot /(norm_c*norm_d));
       }
       else
       {
              TangTemp      =  0;
       }
         
       // Bangle is the bending angle between Tangle is the twisting angle.
       for (int i=0;i<3;i++)  bangle(i+1)  =  BangTemp*result1(i+1);
       for (int i=0;i<3;i++)  tangle(i+1)  =  TangTemp*result2(i+1);
      // cout<<" what is ==  "<< bangle(1)<<bangle(2) << bangle(3)<<endl; 
       //sleep(5000);
       //inertial, bending stiffness, twisting stifness length of each section.
       //HardCode, currently, just set all fibre segments have the same properties.
       Bendinertia    =   PI/8.0*(pow(radiusA,4) + pow(radiusB,4));
       Kb             =   Ey * Bendinertia / L;
       Twistinertia   =   PI * pow(radiusA,4) / 2.0;
       Kt             =   G  * Twistinertia / L;
       //bending and twisting torque
       for(int n = 0; n < 3; n++)
       {
	  btorque(n+1) = Kb*bangle(n+1);
	  ttorque(n+1) = Kt*tangle(n+1);
	  //cout<<" inside function "<< "...torque="<< btorque(n+1) <<endl;	
       }	
       
}  
