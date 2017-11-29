/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "math_extra_liggghts_superquadric.h"
#include "math_extra_liggghts_nonspherical.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#include <cmath>
#include <algorithm>
#if __STDCPP_MATH_SPEC_FUNCS__ >= 201003L
    #define BETA_NAMESPACE std
#elif defined(HAVE_TR1_CMATH)
    #include <tr1/cmath>
    #define BETA_NAMESPACE std::tr1
#elif defined(BOOST_INCLUDED)
    #include "boost/math/special_functions/beta.hpp"
    #define BETA_NAMESPACE boost::math
#else
    
    #include "boost/math/special_functions/beta.hpp"
    #define BETA_NAMESPACE boost::math
    
#endif

#include <map>
#include "math_const.h"

#define PHI_INV 0.61803398874989479

namespace MathExtraLiggghtsNonspherical {

double J4[16];

double check_inequalities(const int i, const double *delta, const double *a, const double *b,
                        const double *A_delta, const double *B_delta, const double *C) {
  double R = 0.0, R0 = 0.0, R1 = 0.0;
  switch (i) {
  case 0: //check A0
    R0 = a[0];
    R1 = b[0]*fabs(C[0]) + b[1]*fabs(C[1])+b[2]*fabs(C[2]);
    R = fabs(A_delta[0]);
    break;
  case 1: //check A1
    R0 = a[1];
    R1 = b[0]*fabs(C[3]) + b[1]*fabs(C[4])+b[2]*fabs(C[5]);
    R = fabs(A_delta[1]);
    break;
  case 2: //check A2
    R0 = a[2];
    R1 = b[0]*fabs(C[6]) + b[1]*fabs(C[7])+b[2]*fabs(C[8]);
    R = fabs(A_delta[2]);
    break;
  case 3: //check B0
    R0 = a[0]*fabs(C[0]) + a[1]*fabs(C[3]) + a[2]*fabs(C[6]);
    R1 = b[0];
    R = fabs(B_delta[0]);
    break;
  case 4: //check B1
    R0 = a[0]*fabs(C[1]) + a[1]*fabs(C[4]) + a[2]*fabs(C[7]);
    R1 = b[1];
    R = fabs(B_delta[1]);
    break;
  case 5: //check B2
    R0 = a[0]*fabs(C[2]) + a[1]*fabs(C[5]) + a[2]*fabs(C[8]);
    R1 = b[2];
    R = fabs(B_delta[2]);
    break;
  case 6: //check //A0xB0
    R0 = a[1]*fabs(C[6]) + a[2]*fabs(C[3]);
    R1 = b[1]*fabs(C[2]) + b[2]*fabs(C[1]);
    R = fabs(C[3]*A_delta[2] - C[6]*A_delta[1]);
    break;
  case 7: //check A0xB1
    R0 = a[1]*fabs(C[7]) + a[2]*fabs(C[4]);
    R1 = b[0]*fabs(C[2]) + b[2]*fabs(C[0]);
    R = fabs(C[4]*A_delta[2] - C[7]*A_delta[1]);
    break;
  case 8: //check A0xB2
    R0 = a[1]*fabs(C[8]) + a[2]*fabs(C[5]);
    R1 = b[0]*fabs(C[1]) + b[1]*fabs(C[0]);
    R = fabs(C[5]*A_delta[2] - C[8]*A_delta[1]);
    break;
  case 9: //check A1xB0
    R0 = a[0]*fabs(C[6]) + a[2]*fabs(C[0]);
    R1 = b[1]*fabs(C[5]) + b[2]*fabs(C[4]);
    R = fabs(C[6]*A_delta[0] - C[0]*A_delta[2]);
    break;
  case 10: //A1xB1
    R0 = a[0]*fabs(C[7]) + a[2]*fabs(C[1]);
    R1 = b[0]*fabs(C[5]) + b[2]*fabs(C[3]);
    R = fabs(C[7]*A_delta[0] - C[1]*A_delta[2]);
    break;
  case 11: //A1xB2
    R0 = a[0]*fabs(C[8]) + a[2]*fabs(C[2]);
    R1 = b[0]*fabs(C[4]) + b[1]*fabs(C[3]);
    R = fabs(C[8]*A_delta[0] - C[2]*A_delta[2]);
    break;
  case 12: //A2xB0
    R0 = a[0]*fabs(C[3]) + a[1]*fabs(C[0]);
    R1 = b[1]*fabs(C[8]) + b[2]*fabs(C[7]);
    R = fabs(C[0]*A_delta[1] - C[3]*A_delta[0]);
    break;
  case 13: //A2xB1
    R0 = a[0]*fabs(C[4]) + a[1]*fabs(C[1]);
    R1 = b[0]*fabs(C[8]) + b[2]*fabs(C[6]);
    R = fabs(C[1]*A_delta[1] - C[4]*A_delta[0]);
    break;
  case 14: //A2xB2
    R0 = a[0]*fabs(C[5]) + a[1]*fabs(C[2]);
    R1 = b[0]*fabs(C[7]) + b[1]*fabs(C[6]);
    R = fabs(C[2]*A_delta[1] - C[5]*A_delta[0]);
    break;
  }
  if(fabs(R0+R1) < 1e-14 and fabs(R) < 1e-14)
    return false;
  else
    return (R > R0+R1);
}

bool obb_intersect(Superquadric *particleA, Superquadric *particleB)
{
  unsigned int inequality_to_start = 0;
  return obb_intersect(particleA, particleB, inequality_to_start);
}

//checking whether bounding boxes intersect or not
// Source: http://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
bool obb_intersect(Superquadric *particleA, Superquadric *particleB, unsigned int &inequality_to_start)
{
  double delta[3];
  LAMMPS_NS::vectorSubtract3D(particleA->center, particleB->center, delta);
  double shape_A[3], shape_B[3];

  for(int k = 0; k < 3; k++) {
    shape_A[k] = particleA->shape[k];
    shape_B[k] = particleB->shape[k];
    particleA->shape[k] *= scale_koef;
    particleB->shape[k] *= scale_koef;
  }

  double C[9];
  MathExtraLiggghtsNonspherical::transpose_times3(particleA->rotation_matrix, particleB->rotation_matrix, C); //C = A^T * B
//check max. 15 inequalities
  double A_delta[3], B_delta[3];
  MathExtraLiggghtsNonspherical::transpose_matvec(particleA->rotation_matrix, delta, A_delta);
  MathExtraLiggghtsNonspherical::transpose_matvec(particleB->rotation_matrix, delta, B_delta);
  bool obb_intersect = true;
  for(unsigned int i = 0; i < 15; i++) {
    int inequality_to_check;
    if(i == 0)
      inequality_to_check = inequality_to_start;
    else if(i == inequality_to_start)
      inequality_to_check = 0;
    else
      inequality_to_check = i;
    if(check_inequalities(inequality_to_check, delta, particleA->shape, particleB->shape, A_delta, B_delta, C)) {
      inequality_to_start = inequality_to_check;
      obb_intersect = false;
      break;
    }
  }
  LAMMPS_NS::vectorCopy3D(shape_A, particleA->shape);
  LAMMPS_NS::vectorCopy3D(shape_B, particleB->shape);
  return obb_intersect;
}

//minimal bounding sphere radius calculation
void bounding_sphere_radius_superquadric(const double *shape, const double *blockiness, double *radius)
{
  if(MathExtraLiggghts::compDouble(blockiness[0],2.0) and MathExtraLiggghts::compDouble(blockiness[1],2.0))
    *radius = std::max(std::max(shape[0], shape[1]), shape[2]);
  else {
    const double a = std::max(shape[0], shape[1]);
    const double b = std::min(shape[0], shape[1]);
    const double c = shape[2];
    if(MathExtraLiggghts::compDouble(blockiness[1],2.0))
      *radius = sqrt(a*a+c*c);
    else {
      const double n1 = std::max(2.1, blockiness[0]);
      const double n2 = std::max(2.1, blockiness[1]);
      double alpha;
      if(!MathExtraLiggghts::compDouble(blockiness[1],2.0))
        alpha = MathExtraLiggghtsNonspherical::pow_abs(b/a, 2.0 / (n2 - 2.0));
      else
        alpha = 0.0;
      const double gamma = pow_abs(1.0 + pow_abs(alpha, n2), n1/n2 - 1.0);
      const double beta = pow_abs(c * c / (a * a) * gamma, 1.0 / (n1 - 2.0));
      const double x0 = 1.0 / pow_abs(pow_abs(1.0 + pow_abs(alpha, n2), n1/n2) + pow_abs(beta, n1), 1.0 / n1);
      const double xtmp = x0 * a;
      const double ytmp = alpha * x0 * b;
      const double ztmp =  beta * x0 * c;
      *radius = sqrt(xtmp*xtmp + ytmp*ytmp + ztmp*ztmp);
    }
  }
  *radius *= scale_koef;
  //*radius = LAMMPS_NS::vectorMag3D(shape);
  /*int index;

  double u_[] = {1.0, 0.0, 0.0};
  for(int i = 0; i <= 90; i = i + 1) {
    double angle = static_cast<double>(i)/180*M_PI;
    double quat_[] = {cos(0.5*angle), 0.0, sin(0.5*angle), 0.0 };
    double shape_[] = {0.354561735, 0.354561735, 0.709123471};
    double aa = shape_[0];
    double bb = shape_[1];
    double cc = shape_[2];
    double blockiness_[] = {6.0, 2.0};
    double p = 1.61;
    double s_ellipsoid = 4.0*M_PI*pow((pow(aa*bb,p) + pow(bb*cc,p) + pow(cc*aa,p))/3.0, 1.0/p);
    double s_cylinder = 2.0 * M_PI*aa*bb + 4.0*M_PI*aa*cc;
    double s_parallelepiped = 8.0*(aa*bb + bb*cc + cc*aa);

    double s_particle;
    area_superquadric(shape_, blockiness_, &s_particle);
    double vol;
    volume_superquadric(shape_, blockiness_, &vol);
    double r = cbrt(0.75/M_PI*vol);
    double s_particle_ortho = crossSectionalArea(u_, shape_, blockiness_, quat_, index);
    double s_sphere = 4.0*M_PI*r*r;
    double s_sphere_ortho = M_PI*r*r;
    double psi = s_sphere / s_particle;
    double psi_ortho = s_sphere_ortho / s_particle_ortho;
    //printf("%3.16e\n", s_particle_ortho);
    printf("%3.16e\n", psi_ortho);
  }
  printf("!");*/
}

// tensor of inertia in particle based coordinate system (diagonal)
//source: http://lrv.fri.uni-lj.si/~franc/SRSbook/geometry.pdf
void inertia_superquadric(const double *shape, const double *blockiness, const double dens, double *ans)
{
  const double eps1 = 2.0 / blockiness[0];
  const double eps2 = 2.0 / blockiness[1];

  const double a1 = shape[0];
  const double a2 = shape[1];
  const double a3 = shape[2];
  const double I_xx = 0.5*a1*a2*a3*eps1*eps2*(a2*a2*BETA_NAMESPACE::beta(1.5*eps2, 0.5*eps2)*BETA_NAMESPACE::beta(0.5*eps1, 2.0*eps1+1.0)+
      4.0*a3*a3*BETA_NAMESPACE::beta(0.5*eps2, 0.5*eps2+1.0)*BETA_NAMESPACE::beta(1.5*eps1, eps1+1.0)) * dens;
  const double I_yy = 0.5*a1*a2*a3*eps1*eps2*(a1*a1*BETA_NAMESPACE::beta(1.5*eps2, 0.5*eps2)*BETA_NAMESPACE::beta(0.5*eps1, 2.0*eps1+1.0)+
      4.0*a3*a3*BETA_NAMESPACE::beta(0.5*eps2, 0.5*eps2+1.0)*BETA_NAMESPACE::beta(1.5*eps1, eps1+1.0)) * dens;
  const double I_zz = 0.5*a1*a2*a3*eps1*eps2*(a1*a1 + a2*a2)*
      BETA_NAMESPACE::beta(1.5*eps2, 0.5*eps2)*BETA_NAMESPACE::beta(0.5*eps1, 2.0*eps1+1.0) * dens;
  ans[0] = I_xx;
  ans[1] = I_yy;
  ans[2] = I_zz;
}

//volume of superquadric particle
//source: http://lrv.fri.uni-lj.si/~franc/SRSbook/geometry.pdf
void volume_superquadric(const double *shape, const double *blockiness, double *result)
{
  const double eps1 = 2.0 / blockiness[0];
  const double eps2 = 2.0 / blockiness[1];
  *result =2.0*shape[0]*shape[1]*shape[2]*eps1*eps2*
      BETA_NAMESPACE::beta(0.5*eps1, eps1 + 1.0)*
      BETA_NAMESPACE::beta(0.5*eps2, 0.5*eps2 + 1.0);
}

void area_superquadric(const double *shape, const double *blockiness, double *result)
{
  double n1 = blockiness[0];
  double n2 = blockiness[1];
  double aa = shape[0];
  double bb = shape[1];
  double cc = shape[2];
  *result = 0.0;

  const int N = 10;
  double sin_angle[N + 1];
  double cos_angle[N + 1];
  double dangle = 0.5*M_PI / static_cast<double>(N); //angle step

  for(int in = 0; in <= N; in++) {
    double angle;
    if(in < N)
      angle = dangle*static_cast<double>(in);
    else
      angle = dangle*(static_cast<double>(in) - 0.01);
    sin_angle[in] = sin(angle); //precompute sinus
    cos_angle[in] = cos(angle); //precompute cosinus
  }

  for(int itheta = 0; itheta < N; itheta++) {
    double sin_theta1 = sin_angle[itheta];
    double cos_theta1 = cos_angle[itheta];
    double sin_theta2 = sin_angle[itheta+1];
    double cos_theta2 = cos_angle[itheta+1];
    double z1 = fabs(sin_theta1);
    double z2 = fabs(sin_theta2);
    for(int iphi = 0; iphi < N; iphi++) {
      double sin_phi1 = sin_angle[iphi];
      double cos_phi1 = cos_angle[iphi];
      double p11[3];
      p11[0] = fabs(cos_phi1)*fabs(cos_theta1);
      p11[1] = fabs(sin_phi1)*fabs(cos_theta1);
      p11[2] = z1;

      double p12[3];
      p12[0] = fabs(cos_phi1)*fabs(cos_theta2);
      p12[1] = fabs(sin_phi1)*fabs(cos_theta2);
      p12[2] = z2;

      double sin_phi2 = sin_angle[iphi+1];
      double cos_phi2 = cos_angle[iphi+1];
      double p21[3];
      p21[0] = fabs(cos_phi2)*fabs(cos_theta1);
      p21[1] = fabs(sin_phi2)*fabs(cos_theta1);
      p21[2] = z1;

      double p22[3];
      p22[0] = fabs(cos_phi2)*fabs(cos_theta2);
      p22[1] = fabs(sin_phi2)*fabs(cos_theta2);
      p22[2] = z2;

      double center_[] = {0,0,0};
      double quat_[] = {1,0,0,0};
      double shape_[] = {aa, bb, cc};
      double blockiness_[] = {n1, n2};
      Superquadric particle(center_, quat_, shape_, blockiness_);
      double f, alpha;

      f = particle.shape_function_local(p11);
      alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/particle.koef + 1.0, 1.0/blockiness[0]);
      for(int i = 0; i < 3; i++)
        p11[i] *= alpha;

      f = particle.shape_function_local(p12);
      alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/particle.koef + 1.0, 1.0/blockiness[0]);
      for(int i = 0; i < 3; i++)
        p12[i] *= alpha;

      f = particle.shape_function_local(p21);
      alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/particle.koef + 1.0, 1.0/blockiness[0]);
      for(int i = 0; i < 3; i++)
        p21[i] *= alpha;

      f = particle.shape_function_local(p22);
      alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/particle.koef + 1.0, 1.0/blockiness[0]);
      for(int i = 0; i < 3; i++)
        p22[i] *= alpha;

      double e1[3], e2[3], e3[4];
      LAMMPS_NS::vectorSubtract3D(p22, p11, e1);
      LAMMPS_NS::vectorSubtract3D(p12, p11, e2);
      LAMMPS_NS::vectorSubtract3D(p21, p11, e3);

      double s1 = LAMMPS_NS::vectorCrossMag3D(e1, e2);
      double s2 = LAMMPS_NS::vectorCrossMag3D(e1, e3);

      *result += (s1 + s2)*4.0;
    }
  }
  /*double result_check1 = 4.0*M_PI*aa*bb; //sphere surface area
  double result_check2 = 8*(aa*bb + bb*cc + cc*aa); //box surface area
  double ecc = sqrt(1.0 - std::min(aa*aa,bb*bb)/std::max(aa*aa,bb*bb));
  double len = 4.0*std::max(aa,bb)*boost::math::ellint_2(ecc);
  double result_check3 = 2.0*M_PI*aa*bb + len*2.0*cc; //cylinder surface area
  printf("%e %e \n", *result, result_check1);*/
}

double crossSectionalArea(const double* u_, const double* shape, const double *blockiness, const double* quat, int& index)
{
  double eps1 = 2.0 / blockiness[0];
  double eps2 = 2.0 / blockiness[1];
  double rot[9];
  MathExtraLiggghtsNonspherical::quat_to_mat(quat, rot);
  double u[3];
  MathExtraLiggghtsNonspherical::transpose_matvec(rot, u_, u);
  LAMMPS_NS::vectorNormalize3D(u);

  if(fabs(u[2]) < 1e-8) {
    u[2] = 0.0;
    LAMMPS_NS::vectorNormalize3D(u);
  }

  double a = shape[0];
  double b = shape[1];
  double c = shape[2];

  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  double result = 0.0;
  int N = 2000;

  if(fabs(uz) > 1e-12) {
    double dphi = M_PI / static_cast<double>(N);
    for(int iphi = 0; iphi < N; iphi ++) {
      double phi =      static_cast<double>(iphi)*dphi;
      double phi_next = static_cast<double>(iphi+1)*dphi;

      double cos_phi = cos(phi);
      double sin_phi = sin(phi);

      double cos_phi_next = cos(phi_next);
      double sin_phi_next = sin(phi_next);

      double tg_theta =      -(ux * pow_abs(cos_phi, 2.0 - eps2) * sign(cos_phi) * (c/a) +
                               uy * pow_abs(sin_phi, 2.0 - eps2) * sign(sin_phi) * (c/b)) / uz;
      double tg_theta_next = -(ux * pow_abs(cos_phi_next, 2.0 - eps2) * sign(cos_phi_next) * (c/a) +
                               uy * pow_abs(sin_phi_next, 2.0 - eps2) * sign(sin_phi_next) * (c/b)) / uz;

      tg_theta = pow_abs(tg_theta, 1.0 / (2.0 - eps1))*sign(tg_theta);
      tg_theta_next = pow_abs(tg_theta_next, 1.0 / (2.0 - eps1))*sign(tg_theta_next);

      double cos_theta2 = 1.0 / (1.0 + tg_theta*tg_theta);
      double sin_theta2 = 1.0 - cos_theta2;
      double cos_theta = sqrt(cos_theta2);
      double sin_theta = sqrt(sin_theta2) * sign(tg_theta);

      double cos_theta_next2 = 1.0 / (1.0 + tg_theta_next*tg_theta_next);
      double sin_theta_next2 = 1.0 - cos_theta_next2;
      double cos_theta_next = sqrt(cos_theta_next2);
      double sin_theta_next = sqrt(sin_theta_next2) * sign(tg_theta_next);

      double p[] =      { a*pow_abs(cos_theta, eps1) * pow_abs(cos_phi, eps2) * sign(cos_theta) * sign(cos_phi),
                          b*pow_abs(cos_theta, eps1) * pow_abs(sin_phi, eps2) * sign(cos_theta) * sign(sin_phi),
                          c*pow_abs(sin_theta, eps1) * sign(sin_theta)};
      double p_next[] = { a*pow_abs(cos_theta_next, eps1) * pow_abs(cos_phi_next, eps2) * sign(cos_theta_next) * sign(cos_phi_next),
                          b*pow_abs(cos_theta_next, eps1) * pow_abs(sin_phi_next, eps2) * sign(cos_theta_next) * sign(sin_phi_next),
                          c*pow_abs(sin_theta_next, eps1) * sign(sin_theta_next)};

      double s[3];
      LAMMPS_NS::vectorCross3D(p, p_next,s);
      double s_ = fabs(LAMMPS_NS::vectorDot3D(s,u));
      result += s_;
    }
    index = 0;
  } else {
    double cos_phi, sin_phi;
    if(fabs(uy) > fabs(ux)) {
      double tg_phi = -ux / uy * b / a;
      tg_phi = pow_abs(tg_phi, 2.0 - eps2) * sign(tg_phi);
      double cos_phi2 = 1.0 / (1.0 + tg_phi*tg_phi);
      cos_phi = sqrt(cos_phi2);
      sin_phi = sqrt(1.0 - cos_phi2)*sign(tg_phi);
      index = 1;
    } else {
      double ctg_phi = -uy / ux * a / b;
      ctg_phi = pow_abs(ctg_phi, 2.0 - eps2) * sign(ctg_phi);
      double sin_phi2 = 1.0 / (1.0 + ctg_phi*ctg_phi);
      sin_phi = sqrt(sin_phi2);
      cos_phi = sqrt(1.0 - sin_phi2)*sign(ctg_phi);
      index = 2;
    }
    double dtheta = M_PI / static_cast<double>(N);
    for(int itheta = 0; itheta < N; itheta++) {
      double theta =      dtheta*static_cast<double>(itheta);
      double theta_next = dtheta*static_cast<double>(itheta+1);
      double cos_theta = cos(theta);
      double sin_theta = sin(theta);

      double cos_theta_next = cos(theta_next);
      double sin_theta_next = sin(theta_next);

      double p[] =      { a * pow_abs(sin_theta, eps1) * pow_abs(cos_phi,eps2) * sign(sin_theta) * sign(cos_phi),
                          b * pow_abs(sin_theta, eps1) * pow_abs(sin_phi,eps2) * sign(sin_theta) * sign(sin_phi),
                          c * pow_abs(cos_theta, eps1) * sign(cos_theta)};
      double p_next[] = { a * pow_abs(sin_theta_next, eps1) * pow_abs(cos_phi,eps2) * sign(sin_theta_next) * sign(cos_phi),
                          b * pow_abs(sin_theta_next, eps1) * pow_abs(sin_phi,eps2) * sign(sin_theta_next) * sign(sin_phi),
                          c * pow_abs(cos_theta_next, eps1) * sign(cos_theta_next)};

      double s[3];
      LAMMPS_NS::vectorCross3D(p, p_next,s);
      double s_ = fabs(LAMMPS_NS::vectorDot3D(s,u));
      result += s_;
    }
  }
  return result;
}

void get_point(int iphi, int nphi, int itheta, int ntheta, Superquadric &particle, std::map<std::pair<int,int>, point> &point_map, double *p) {
  std::pair<int, int> pair;
  pair = std::make_pair(iphi, itheta);
  std::map<std::pair<int,int>, point>::iterator it = point_map.find(pair);
  if(it == point_map.end()) {
    particle.reference_point(iphi, nphi, itheta, ntheta, p);
    point_map[pair] = point(p);
  } else {
    point cur_map_element = it->second;
    p[0] = cur_map_element.x;
    p[1] = cur_map_element.y;
    p[2] = cur_map_element.z;
  }
}

/* OBSOLETE CODE
TODO: remove?
void initial_estimate(SurfacesIntersectData &stdata, int *nphi_start, int *ntheta_start,
    const double *point_i_estimate, const double *point_j_estimate, double *result_pointA, double *result_pointB)
{
  #ifdef SUPERQUADRIC_ACTIVE_FLAG

  std::map<std::pair<int,int>, point> map1;
  std::map<std::pair<int,int>, point> map2;

  double p1[3], p2[3];
  double p1_[3], p2_[3];
  int nphi = *nphi_start;
  int ntheta = *ntheta_start;
  double distance_first_estimate = LAMMPS_NS::pointDistance(point_i_estimate, point_j_estimate);

  int iphi1, itheta1, iphi2, itheta2;
  Superquadric particleA(stdata.pos_i, stdata.quat_i, stdata.shape_i, stdata.blockiness_i);
  Superquadric particleB(stdata.pos_j, stdata.quat_j, stdata.shape_j, stdata.blockiness_j);

  particleA.pre_initial_estimate(point_i_estimate, nphi, ntheta, &iphi1, &itheta1);
  particleB.pre_initial_estimate(point_j_estimate, nphi, ntheta, &iphi2, &itheta2);

  get_point(iphi1, nphi, itheta1, ntheta, particleA, map1, p1);
  get_point(iphi2, nphi, itheta2, ntheta, particleB, map2, p2);

  bool next_iter = true;
  int it = 0;
  while(next_iter) {
    next_iter = false;
    get_point(iphi1, nphi, itheta1, ntheta, particleA, map1, p1);
    get_point(iphi2, nphi, itheta2, ntheta, particleB, map2, p2);
    it ++;
    int ir = 0;
    int shifts[8][2] = {{ 0,  1},
                        { 0, -1},
                        { 1,  0},
                        {-1,  0},
                        { 1,  1},
                        { 1, -1},
                        {-1,  1},
                        {-1, -1}};
    double rsq = LAMMPS_NS::pointDistance(p1, p2);
    int iphi1_ = iphi1;
    int itheta1_ = itheta1;
    int iphi2_ = iphi2;
    int itheta2_ = itheta2;
    for(int ishift1 = 0; ishift1 < 8; ishift1++) {
      double p1_new[3], p2_new[3];
      int iphi1_new =     (iphi1 + shifts[ishift1][0] +   nphi - 1) % (nphi - 1);
      int itheta1_new = (itheta1 + shifts[ishift1][1] + ntheta - 1) % (ntheta - 1);
      get_point(iphi1_new, nphi, itheta1_new, ntheta, particleA, map1, p1_new);
      for(int ishift2 = 0; ishift2 < 8; ishift2++) {
        ir ++;
        int iphi2_new =     (iphi2 + shifts[ishift2][0] +   nphi - 1) % (nphi - 1);
        int itheta2_new = (itheta2 + shifts[ishift2][1] + ntheta - 1) % (ntheta - 1);
        get_point(iphi2_new, nphi, itheta2_new, ntheta, particleB, map2, p2_new);
        double rsq_ = LAMMPS_NS::pointDistance(p1_new, p2_new);
        if(rsq_ < rsq) {
          rsq = rsq_;
          iphi1_ = iphi1_new;
          itheta1_ = itheta1_new;
          iphi2_ = iphi2_new;
          itheta2_ = itheta2_new;
          LAMMPS_NS::vectorCopy3D(p1_new, p1_);
          LAMMPS_NS::vectorCopy3D(p2_new, p2_);
          next_iter = true;
        }
      }
    }
    iphi1 = iphi1_;
    itheta1 = itheta1_;
    iphi2 = iphi2_;
    itheta2 = itheta2_;
    LAMMPS_NS::vectorCopy3D(p1_, p1);
    LAMMPS_NS::vectorCopy3D(p2_, p2);
  }

  double rsq = LAMMPS_NS::pointDistance(p1, p2);
  //double r = sqrt(rsq);
  if(distance_first_estimate < rsq) {
    LAMMPS_NS::vectorCopy3D(point_i_estimate, result_pointA);
    LAMMPS_NS::vectorCopy3D(point_j_estimate, result_pointB);
  } else {
    LAMMPS_NS::vectorCopy3D(p1, result_pointA);
    LAMMPS_NS::vectorCopy3D(p2, result_pointB);
  }

  #endif
}*/

double calc_F(Superquadric *particleA, Superquadric *particleB, double &f1, double &f2, double *gradA, double *gradB, double *hess1, double *hess2, const double *point, double mu, double *F, double *merit)
{
  /*particleA->shape_function_gradient_global(point, gradA);
  particleB->shape_function_gradient_global(point, gradB);
  f1 = particleA->shape_function_global(point);
  f2 = particleB->shape_function_global(point);*/

  particleA->shape_function_props_global(point, &f1, gradA, hess1);
  particleB->shape_function_props_global(point, &f2, gradB, hess2);

  double mu_sq = mu*mu;
  F[0] = gradA[0] + mu_sq*gradB[0];
  F[1] = gradA[1] + mu_sq*gradB[1];
  F[2] = gradA[2] + mu_sq*gradB[2];
  F[3] = f1 - f2;

  double v[3];
  MathExtra::cross3(gradA, gradB, v);
  double sine = (fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * LAMMPS_NS::vectorDot3D(gradB, gradB))));
  double aaa = fabs(f1 - f2)/(std::max(fabs(f1),fabs(f2))+1e-16);
  *merit = std::max(sine, aaa*aaa);
  return F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
}

double calc_F(Superquadric *particleA, Superquadric *particleB, double &f1, double &f2, double *gradA, double *gradB, double *hess1, double *hess2, const double *point, double *mu, double *F, double *merit)
{
  /*particleA->shape_function_gradient_global(point, gradA);
  particleB->shape_function_gradient_global(point, gradB);
  f1 = particleA->shape_function_global(point);
  f2 = particleB->shape_function_global(point);*/

  particleA->shape_function_props_global(point, &f1, gradA, hess1);
  particleB->shape_function_props_global(point, &f2, gradB, hess2);

  double mu_sq = -LAMMPS_NS::vectorDot3D(gradA,gradB) / LAMMPS_NS::vectorDot3D(gradB,gradB);
  *mu = sqrt(fabs(mu_sq));
  F[0] = gradA[0] + mu_sq*gradB[0];
  F[1] = gradA[1] + mu_sq*gradB[1];
  F[2] = gradA[2] + mu_sq*gradB[2];
  F[3] = f1 - f2;

  double v[3];
  MathExtra::cross3(gradA, gradB, v);
  double sine = (fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * LAMMPS_NS::vectorDot3D(gradB, gradB))));
  double aaa = fabs(f1 - f2)/(std::max(fabs(f1),fabs(f2))+1e-16);
  *merit = std::max(sine, aaa*aaa);
  return F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
}

double calc_F_new(Superquadric *particle1, Superquadric *particle2, const double *x2, double *gradA, double *gradB, const double alpha, const double beta, double *F, double *merit)
{
  double x1[3];
  double f1, f2;

  particle2->shape_function_props_global(x2, &f2, gradB, NULL);
  x1[0] = x2[0] + beta*gradB[0];
  x1[1] = x2[1] + beta*gradB[1];
  x1[2] = x2[2] + beta*gradB[2];
  particle1->shape_function_props_global(x1, &f1, gradA, NULL);

  F[0] = gradA[0]+alpha*gradB[0];
  F[1] = gradA[1]+alpha*gradB[1];
  F[2] = gradA[2]+alpha*gradB[2];
  F[3] = f1;
  F[4] = f2;

  double v[3];
  MathExtra::cross3(gradA, gradB, v);
  double sin_ = sqrt(fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * LAMMPS_NS::vectorDot3D(gradB, gradB))));
  *merit = fabs(sin_) + fabs(f1) + fabs(f2);

  return sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3] + F[4]*F[4]);
}

double calc_F_new(Superquadric *particle1, Superquadric *particle2, double *x1, const double *x2, double *gradA, double *gradB, const double alpha, const double beta, double *F, double *merit)
{
  double f1, f2;
  particle2->shape_function_props_global(x2, &f2, gradB, NULL);
  x1[0] = x2[0] + beta*gradB[0];
  x1[1] = x2[1] + beta*gradB[1];
  x1[2] = x2[2] + beta*gradB[2];
  particle1->shape_function_props_global(x1, &f1, gradA, NULL);

  F[0] = gradA[0]+alpha*gradB[0];
  F[1] = gradA[1]+alpha*gradB[1];
  F[2] = gradA[2]+alpha*gradB[2];
  F[3] = f1;
  F[4] = f2;

  double v[3];
  MathExtra::cross3(gradA, gradB, v);
  double sin_ = sqrt(fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * LAMMPS_NS::vectorDot3D(gradB, gradB))));
  *merit = fabs(sin_) + fabs(f1) + fabs(f2);

  return sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3] + F[4]*F[4]);
}

double invf(int i,int j,const double* m){

    int o = 2+(j-i);

    i += 4+o;
    j += 4-o;

    #define e(a,b) m[ ((j+b)%4)*4 + ((i+a)%4) ]

    double inv =
     + e(+1,-1)*e(+0,+0)*e(-1,+1)
     + e(+1,+1)*e(+0,-1)*e(-1,+0)
     + e(-1,-1)*e(+1,+0)*e(+0,+1)
     - e(-1,-1)*e(+0,+0)*e(+1,+1)
     - e(-1,+1)*e(+0,-1)*e(+1,+0)
     - e(+1,-1)*e(-1,+0)*e(+0,+1);

    return (o%2)?inv : -inv;

    #undef e
}

double inverseMatrix4x4(const double *m, double *out)
{

    double inv[16];

    for(int i = 0;i < 4; i++)
        for(int j = 0; j < 4; j++)
            inv[j*4+i] = invf(i,j,m);

    double D = 0;

    for(int k = 0; k < 4; k++) D += m[k] * inv[k*4];

    if (D == 0) return D;

    double D_inv = 1.0 / D;

    for (int i = 0; i < 16; i++)
        out[i] = inv[i] * D_inv;

    return D;
}

double determinant_4x4(double *mat)
{
  const double m00 = mat[0]; const double m01 = mat[1]; const double m02 = mat[2]; const double m03 = mat[3];
  const double m10 = mat[4]; const double m11 = mat[5]; const double m12 = mat[6]; const double m13 = mat[7];
  const double m20 = mat[8]; const double m21 = mat[9]; const double m22 = mat[10]; const double m23 = mat[11];
  const double m30 = mat[12]; const double m31 = mat[13]; const double m32 = mat[14]; const double m33 = mat[15];

  double value =
     m03*m12*m21*m30 - m02*m13*m21*m30 - m03*m11*m22*m30 + m01*m13*m22*m30+
     m02*m11*m23*m30 - m01*m12*m23*m30 - m03*m12*m20*m31 + m02*m13*m20*m31+
     m03*m10*m22*m31 - m00*m13*m22*m31 - m02*m10*m23*m31 + m00*m12*m23*m31+
     m03*m11*m20*m32 - m01*m13*m20*m32 - m03*m10*m21*m32 + m00*m13*m21*m32+
     m01*m10*m23*m32 - m00*m11*m23*m32 - m02*m11*m20*m33 + m01*m12*m20*m33+
     m02*m10*m21*m33 - m00*m12*m21*m33 - m01*m10*m22*m33 + m00*m11*m22*m33;
  return value;
}

void calc_contact_point(Superquadric *particleA, Superquadric *particleB,
    double ratio, const double *initial_point1, double *result_point, double &fi, double &fj, bool *fail, LAMMPS_NS::Error *error)
{
  const double tol1 = 1e-10; //tolerance
  const double tol2 = 1e-12;
  *fail = false;

  double mu, mu_sq;
  double F[4], F1[4], F2[4];
  double point[3];

  double mu1, mu2;
  double gradA1[3], gradB1[3];

  double fi1, fj1, fi2, fj2;

  double merit01, merit02, merit0;

  double res01 = calc_F(particleA, particleB, fi1, fj1, gradA1, gradB1, NULL, NULL, initial_point1, &mu1, F1, &merit01);
  if(merit01 < tol1) {
    LAMMPS_NS::vectorCopy3D(initial_point1, result_point);
    LAMMPS_NS::vectorCopy3D(gradA1, particleA->gradient);
    LAMMPS_NS::vectorCopy3D(gradB1, particleB->gradient);
    fi = fi1;
    fj = fj1;
    return; //finish if initial guess is good enough
  }

  double initial_point2[3];
  for(int k = 0; k < 3; k++)
    initial_point2[k] = ratio*particleB->center[k] + (1.0 - ratio)*particleA->center[k]; //solution for spheres
  double gradA2[3], gradB2[3];
  double res02 = calc_F(particleA, particleB, fi2, fj2, gradA2, gradB2, NULL, NULL, initial_point2, &mu2, F2, &merit02);
  double res0;
  if(res01 < res02) {
    res0 = res01; //solution from previous step is better than for spheres
    LAMMPS_NS::vectorCopy3D(gradA1, particleA->gradient);
    LAMMPS_NS::vectorCopy3D(gradB1, particleB->gradient);
    LAMMPS_NS::vectorCopy3D(initial_point1, point);
    LAMMPS_NS::vectorCopy4D(F1, F);
    fi = fi1;
    fj = fj1;
    mu = mu1;
    merit0 = merit01;
  } else { //solution for spheres better than from the previous step one
    res0 = res02;
    LAMMPS_NS::vectorCopy3D(gradA2, particleA->gradient);
    LAMMPS_NS::vectorCopy3D(gradB2, particleB->gradient);
    LAMMPS_NS::vectorCopy3D(initial_point2, point);
    LAMMPS_NS::vectorCopy4D(F2, F);
    fi = fi2;
    fj = fj2;
    mu = mu2;
    merit0 = merit02;
  }

  double point0[3];
  LAMMPS_NS::vectorCopy3D(point, point0);

  if(merit0 < tol1) {
    LAMMPS_NS::vectorCopy3D(point, result_point);
    return; //finish if the solution for spheres is good enough
  }

  double merit1 = merit0;
  double merit2 = merit0;

  double res1 = res0;
  double res2 = res1;

  double size_i = MathExtraLiggghts::min(particleA->shape[0],particleA->shape[1],particleA->shape[2]);
  double size_j = MathExtraLiggghts::min(particleB->shape[0],particleB->shape[1],particleB->shape[2]);

  double size = std::min(size_i, size_j);

  J4[15] = 0.0;
  double delta[4], delta_0[4];
  zeros(delta, 4);
  zeros(delta_0, 4);
  const int Niter = 100000;
  double pointb[3], pointa[3];
  double J4_inv[16];

  for(int iter = 0; iter < Niter; iter++) {

    merit2 = merit1;
    res2 = res1;

    particleA->shape_function_hessian_global(point, particleA->hessian);
    particleB->shape_function_hessian_global(point, particleB->hessian);

    mu_sq = mu*mu;
    for(int i = 0; i < 3; i++) { //construct Jacobian
      for(int j = 0; j < 3; j++) {
        J4[4*i + j] = particleA->hessian[i*3+j] + mu_sq*particleB->hessian[i*3+j];
      }
      J4[3+4*i] = 2.0*mu*particleB->gradient[i];
      J4[i+4*3] = particleA->gradient[i] - particleB->gradient[i];
    }

    double det = determinant_4x4(J4);
    if(fabs(det) > 1.0) {
      inverseMatrix4x4(J4, J4_inv);
      MathExtraLiggghtsNonspherical::matvecN(J4_inv, F, 4, delta);
    } else {
      vectorCopyN(delta, 4, delta_0);
      GMRES<4,4>(J4, F, delta_0, delta); //solve linear system
    }

    bool nan_found = false;
    for(unsigned int i = 0; i < 4; i++) {
      if(std::isnan(delta[i])) {
        nan_found = true; //check if there are NaNs in the solution
#ifdef LIGGGHTS_DEBUG
        error->warning(FLERR,"NaN found in calc_contact_point while solving linear system!\n");
        printf_debug_data(particleA, particleB, point0, error);
#endif
        *fail = true;
        break;
      }
    }

    bool converged = false;
    double deltax = MathExtraLiggghts::max(fabs(delta[0]),fabs(delta[1]),fabs(delta[2]));
    if(deltax == 0.0) {
      LAMMPS_NS::vectorCopy3D(point, result_point);
      LAMMPS_NS::vectorCopy3D(point, result_point);
      converged = true;
    } else {
    if(nan_found)
      break;
    else {
      double point_[3], mu_;
      point_[0] = point[0] - delta[0];
      point_[1] = point[1] - delta[1];
      point_[2] = point[2] - delta[2];
      mu_ = mu - delta[3];
      double fi_, fj_;
        double merit2_;
        double res2_ = calc_F(particleA, particleB, fi_, fj_, particleA->gradient, particleB->gradient, NULL, NULL, point_, mu_, F, &merit2_);
        if(res2_ < res1 or merit2_ < tol1 or deltax < tol2 * size) {
          merit2 = merit2_;
          res2 = res2_;
        mu = mu_;
        fi = fi_;
        fj = fj_;
        LAMMPS_NS::vectorCopy3D(point_, point);
          if(merit2 < tol1 or deltax < tol2 * size)
          converged = true;
      } else { //make the steepest descent of the residual
        double a = 0.0;
        double b = 1.0;
        double eps = fabs(b-a);

        double res2a, res2b;
          double merit2a, merit2b;

        while(eps > 1e-8) {
          double alpha1 = b - (b - a)*PHI_INV;
          double alpha2 = a + (b - a)*PHI_INV;

          yabx3D(point, -alpha1, delta, pointa);
          mu1 = mu - alpha1*delta[3];
            res2a = calc_F(particleA, particleB, fi1, fj1, gradA1, gradB1, NULL, NULL, pointa, mu1, F1, &merit2a);

          yabx3D(point, -alpha2, delta, pointb);
          mu2 = mu - alpha2*delta[3];
            res2b = calc_F(particleA, particleB, fi2, fj2, gradA2, gradB2, NULL, NULL, pointb, mu2, F2, &merit2b);
          if(std::min(res2a, res2b) < res1) {
            if(res2a < res2b) {
                res2 = res2a;
              LAMMPS_NS::vectorCopy4D(F1, F);
              mu = mu1;
                LAMMPS_NS::vectorCopy3D(gradA1, particleA->gradient);
                LAMMPS_NS::vectorCopy3D(gradB1, particleB->gradient);
              LAMMPS_NS::vectorCopy3D(pointa, point);
              fi = fi1;
              fj = fj1;
                merit2 = merit2a;
            } else {
                res2 = res2b;
              LAMMPS_NS::vectorCopy4D(F2, F);
              mu = mu2;
                LAMMPS_NS::vectorCopy3D(gradA2, particleA->gradient);
                LAMMPS_NS::vectorCopy3D(gradB2, particleB->gradient);
              LAMMPS_NS::vectorCopy3D(pointb, point);
              fi = fi2;
              fj = fj2;
                merit2 = merit2b;
            }
              if(merit2 < tol1)
              converged = true;
            eps = 0.0;
          } else  {
            if(res2a > res2b)
              a = alpha1;
            else
              b = alpha2;
            eps = fabs(b-a);
          }
        }
        if(eps != 0.0) {
          LAMMPS_NS::vectorCopy4D(F2, F);
          mu = mu2;
            LAMMPS_NS::vectorCopy3D(gradA2, particleA->gradient);
            LAMMPS_NS::vectorCopy3D(gradB2, particleB->gradient);
          LAMMPS_NS::vectorCopy3D(pointb, point);
          fi = fi2;
          fj = fj2;
#ifdef LIGGGHTS_DEBUG
            if(std::min(merit2a,merit2b)>merit1) {
            error->warning(FLERR, "Could not find optimal step by Golden section algorithm, epsilon<1e-6");
              printf_debug_data(particleA, particleB, point0, error);
              printf("res2b: %e merit: %e\n",res2b, merit2b);
          }
#endif
            res2 = res2b;
            merit2 = merit2b;
            if(merit2 < tol1)
            converged = true;
        }
      }
    }

      if(fabs(res1-res2)/std::max(res1, res2) < 1e-3)
        converged = true;
    }

    merit1 = merit2;
    res1 = res2;
    if(iter >= Niter-1) {
        *fail = true;
#ifdef LIGGGHTS_DEBUG
        error->warning(FLERR, "Contact detection algorithm could not converge to desired tolerance in desired number of iterations\n");
        printf_debug_data(particleA, particleB, point0, error);
#endif
    }
    if(converged || *fail) {
      LAMMPS_NS::vectorCopy3D(point, result_point);
      break;
    }
  }
}

const int ndim = 8;

double calc_F_IP(Superquadric *particleA, Superquadric *particleB,
    double *pointA, double *pointB,
    double *gradA, double *gradB,
    double lA, double lB,
    double *F, double *merit)
{
  double fA, fB;

  particleA->shape_function_props_global(pointA, &fA, gradA, NULL);
  particleB->shape_function_props_global(pointB, &fB, gradB, NULL);

  double delta[3];
  LAMMPS_NS::vectorSubtract3D(pointA, pointB, delta);

  F[0] = delta[0] + lA*gradA[0];
  F[1] = delta[1] + lA*gradA[1];
  F[2] = delta[2] + lA*gradA[2];

  F[3] = gradA[0] + lB*gradB[0];
  F[4] = gradA[1] + lB*gradB[1];
  F[5] = gradA[2] + lB*gradB[2];

  F[6] = fA;
  F[7] = fB;

  double vA[3], vB[3];
  MathExtra::cross3(gradA, delta, vA);
  MathExtra::cross3(gradB, delta, vB);
  double sin2A = fabs(LAMMPS_NS::vectorDot3D(vA,vA)) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * (LAMMPS_NS::vectorDot3D(delta, delta)+1e-16));
  double sin2B = fabs(LAMMPS_NS::vectorDot3D(vB,vB)) / (LAMMPS_NS::vectorDot3D(gradB, gradB) * (LAMMPS_NS::vectorDot3D(delta, delta)+1e-16));

  *merit = fabs(sin2A) + fabs(sin2B) + fabs(fA) + fabs(fB);
  double res = 0.0;
  for(int i = 0; i < ndim; i++)
    res += F[i]*F[i];
  return res;
  }

double calc_F_IP(Superquadric *particleA, Superquadric *particleB,
    const double *pointA,  const double *pointB,
    double *gradA, double *gradB,
    double *lA_, double *lB_, double *F, double *merit, double *sine)
{
  double fA, fB;

  particleA->shape_function_props_global(pointA, &fA, gradA, NULL);
  particleB->shape_function_props_global(pointB, &fB, gradB, NULL);

  double delta[3];
  LAMMPS_NS::vectorSubtract3D(pointA, pointB, delta);

  *lA_ = -MathExtra::dot3(gradA, delta) / MathExtra::dot3(gradA, gradA);
  *lB_ = -MathExtra::dot3(gradB, gradA) / MathExtra::dot3(gradB, gradB);

  double lA = *lA_;
  double lB = *lB_;

  F[0] = delta[0] + lA*gradA[0];
  F[1] = delta[1] + lA*gradA[1];
  F[2] = delta[2] + lA*gradA[2];

  F[3] = gradA[0] + lB*gradB[0];
  F[4] = gradA[1] + lB*gradB[1];
  F[5] = gradA[2] + lB*gradB[2];

  F[6] = fA;
  F[7] = fB;

  double vA[3], vB[3];
  MathExtra::cross3(gradA, delta, vA);
  MathExtra::cross3(gradB, delta, vB);

  double sin2A = fabs(LAMMPS_NS::vectorDot3D(vA,vA)) / (LAMMPS_NS::vectorDot3D(gradA, gradA) * (LAMMPS_NS::vectorDot3D(delta, delta)+1e-16));
  double sin2B = fabs(LAMMPS_NS::vectorDot3D(vB,vB)) / (LAMMPS_NS::vectorDot3D(gradB, gradB) * (LAMMPS_NS::vectorDot3D(delta, delta)+1e-16));

  *merit = fabs(sin2A) + fabs(sin2B) + fabs(fA) + fabs(fB);
  *sine = 0.5*sqrt(fabs(sin2A) + fabs(sin2B));

  double res = 0.0;
  for(int i = 0; i < ndim; i++)
    res += F[i]*F[i];
  return res;
    }

//a minimal distance algorithm based on a common normal concept
double minimal_distance(Superquadric *particleA, Superquadric *particleB, const double *initial_pointA_1, const double *initial_pointB_1,
                                                                              double *result_pointA,         double *result_pointB,
                                                                               bool *fail)
{
  *fail = false;

  double tol0 = 1e-26;
  double tol1 = 1e-12;
  double tol2 = 1e-16;
  double pointA[3], pointB[3];

  double F[ndim];
  double lA, lB;

  double F1[ndim], gradA1[3], gradB1[3];
  double lA1, lB1, merit1;

  double F2[ndim], gradA2[3], gradB2[3];
  double lA2, lB2, merit2;

  double merit0;

  double sine1;

 //calculate the residual vector for the first pair candidate
  double res0;
  double merit01;
  //calculate the residual vector for the first pair candidate
  double res01 = calc_F_IP(particleA, particleB, initial_pointA_1, initial_pointB_1, particleA->gradient, particleB->gradient, &lA1, &lB1, F1, &merit01, &sine1);
  if(merit01 < tol1 or res01 < tol0) {
    LAMMPS_NS::vectorCopy3D(initial_pointA_1, result_pointA);
    LAMMPS_NS::vectorCopy3D(initial_pointB_1, result_pointB);
    return merit01;
  }

  vectorCopyN(F1, ndim, F);
  LAMMPS_NS::vectorCopy3D(initial_pointA_1, pointA);
  LAMMPS_NS::vectorCopy3D(initial_pointB_1, pointB);
  lA = lA1;
  lB = lB1;
  res0 = res01;
  merit0 = merit01;

  merit1 = merit0;
  merit2 = merit0;

  double pointA0[3],pointB0[3];
  LAMMPS_NS::vectorCopy3D(pointA, pointA0);
  LAMMPS_NS::vectorCopy3D(pointB, pointB0);

  double delta[ndim], delta_0[ndim];
  zeros(delta, ndim);
  for(int i = 0; i < ndim; i++)
    delta_0[i] = 1.0;
  double hessA[9], hessB[9];

  double res1 = res0;
  double res2 = res1;

  bool nan_found = false;

  bool converged = false;

  int Niter = 10000;

  double sizeA = MathExtraLiggghts::min(particleA->shape[0],particleA->shape[1],particleA->shape[2]);
  double sizeB = MathExtraLiggghts::min(particleB->shape[0],particleB->shape[1],particleB->shape[2]);

  double size = std::min(sizeA, sizeB);

  for(int iter = 0; iter < Niter; iter ++)
  {

    res2 = res1;
    merit2 = merit1;

    particleA->shape_function_hessian_global(pointA, hessA);
    particleB->shape_function_hessian_global(pointB, hessB);

    double Jacobian[] =
    {
        1.0+lA*hessA[0],     lA*hessA[1],    lA*hessA[2],-1.0,  0.0,  0.0,   particleA->gradient[0], 0,
            lA*hessA[3], 1.0+lA*hessA[4],    lA*hessA[5], 0.0, -1.0,  0.0,   particleA->gradient[1], 0,
            lA*hessA[6],     lA*hessA[7],1.0+lA*hessA[8], 0.0,  0.0, -1.0,   particleA->gradient[2], 0,

        hessA[0], hessA[1], hessA[2], lB*hessB[0], lB*hessB[1], lB*hessB[2], 0.0, particleB->gradient[0],
        hessA[3], hessA[4], hessA[5], lB*hessB[3], lB*hessB[4], lB*hessB[5], 0.0, particleB->gradient[1],
        hessA[6], hessA[7], hessA[8], lB*hessB[6], lB*hessB[7], lB*hessB[8], 0.0, particleB->gradient[2],

        particleA->gradient[0], particleA->gradient[1], particleA->gradient[2], 0.0, 0.0, 0.0,  0.0, 0.0,
        0.0, 0.0, 0.0, particleB->gradient[0], particleB->gradient[1], particleB->gradient[2],  0.0, 0.0

    };

    vectorCopyN(delta, ndim, delta_0);
    GMRES<ndim,ndim>(Jacobian, F, delta_0, delta); //solve linear system

    double deltax1 = MathExtraLiggghts::max(fabs(delta[0]),fabs(delta[1]),fabs(delta[2]));
    double deltax2 = MathExtraLiggghts::max(fabs(delta[3]),fabs(delta[4]),fabs(delta[5]));

    double deltax = std::max(deltax1, deltax2);
    if(deltax == 0.0) {
      LAMMPS_NS::vectorCopy3D(pointA, result_pointA);
      LAMMPS_NS::vectorCopy3D(pointB, result_pointB);

      converged = true;
    } else {

      for(int i = 0; i < ndim; i++) {
      if(std::isnan(delta[i])) {
        nan_found = true;
        break;
      }
    }
    if(nan_found)
      break;
      if(deltax > 0.01*size) {
        double mag = normN(delta, ndim);
        for(int i = 0; i < ndim; i++)
          delta[i] *= 0.01*size/mag;
      }

      double lA_, lB_;
      double pointA_[3],pointB_[3];
      double pointA1[3], pointB1[3];
      double pointA2[3], pointB2[3];
      LAMMPS_NS::vectorSubtract3D(pointA, delta, pointA_);
      LAMMPS_NS::vectorSubtract3D(pointB, delta+3, pointB_);
      lA_ = lA - delta[6];
      lB_ = lB - delta[7];

      double merit2_;
      double res2_ = calc_F_IP(particleA, particleB, pointA_, pointB_, particleA->gradient, particleB->gradient, lA_, lB_, F, &merit2_);
      if(res2_ < res1 or (res2_ < tol0 or merit2_ < tol1) or  deltax < tol2 * size) {
        res2 = res2_;
        merit2 = merit2_;
        LAMMPS_NS::vectorCopy3D(pointA_, pointA);
        LAMMPS_NS::vectorCopy3D(pointB_, pointB);

        lA = lA_;
        lB = lB_;

        if((res2_ < tol0 or merit2_ < tol1) || deltax < tol2 * size)
        converged = true;
      } else { //steepest descent
        double a = 0.0;
        double b = 1.0;
        double eps = fabs(b-a);
      double res2a, res2b;
        double merit2a, merit2b;

        while(eps > 1e-12) { //the golden section algorithm
          double alpha1 = b - (b - a)*PHI_INV;
          double alpha2 = a + (b - a)*PHI_INV;
          yabx3D(pointA, -alpha1, delta, pointA1);
          yabx3D(pointB, -alpha1, delta+3, pointB1);
          lA1 = lA - alpha1 * delta[6];
          lB1 = lB - alpha1 * delta[7];

          res2a = calc_F_IP(particleA, particleB, pointA1, pointB1, gradA1, gradB1, lA1, lB1, F1, &merit2a);

          yabx3D(pointA, -alpha2, delta, pointA2);
          yabx3D(pointB, -alpha2, delta+3, pointB2);
          lA2 = lA - alpha2 * delta[6];
          lB2 = lB - alpha2 * delta[7];

          res2b = calc_F_IP(particleA, particleB, pointA2, pointB2, gradA2, gradB2, lA2, lB2, F2, &merit2b);

        if(std::min(res2a, res2b) < res1) {
          if(res2a < res2b) {
              res2 = res2a;
              vectorCopyN(F1, ndim, F);
              LAMMPS_NS::vectorCopy3D(gradA1, particleA->gradient);
              LAMMPS_NS::vectorCopy3D(gradB1, particleB->gradient);
              LAMMPS_NS::vectorCopy3D(pointA1, pointA);
              LAMMPS_NS::vectorCopy3D(pointB1, pointB);
              lA = lA1;
              lB = lB1;

              merit2 = merit2a;
          } else {
              res2 = res2b;
              vectorCopyN(F2, ndim, F);
              LAMMPS_NS::vectorCopy3D(gradA2, particleA->gradient);
              LAMMPS_NS::vectorCopy3D(gradB2, particleB->gradient);
              LAMMPS_NS::vectorCopy3D(pointA2, pointA);
              LAMMPS_NS::vectorCopy3D(pointB2, pointB);
              lA = lA2;
              lB = lB2;

              merit2 = merit2b;
          }
            if(res2b < tol0 or merit1 < tol1)
            converged = true;
          eps = 0.0;
        } else {
          if(res2a > res2b)
              a = alpha1;
          else
              b = alpha2;
            eps = fabs(b-a);
         }
       }
        if(eps != 0.0)
          converged = true;
      }
    }

    if(fabs(res1 - res2)/std::max(res1, res2) < 1e-3)
      converged = true;

    res1 = res2;
    merit1 = merit2;

    if(converged || iter == Niter - 1) {
      LAMMPS_NS::vectorCopy3D(pointA, result_pointA);
      LAMMPS_NS::vectorCopy3D(pointB, result_pointB);
      break;
    }
  }
  return merit1;
}

bool capsules_intersect(Superquadric *particleA, Superquadric *particleB, double *capsule_contact_point) {

  int main_axis_i = 2, main_axis_j = 2; //z-axis is main axis for cylinders
  if(!particleA->isCylinder) {
    if(particleA->shape[0] >= std::max(particleA->shape[1], particleA->shape[2]))
      main_axis_i = 0;
    else if(particleA->shape[1] >= std::max(particleA->shape[0], particleA->shape[2]))
      main_axis_i = 1;
    else
      main_axis_i = 2;
    //printf("!");
  }

  if(!particleB->isCylinder) {
    if(particleB->shape[0] >= std::max(particleB->shape[1], particleB->shape[2]))
      main_axis_j = 0;
    else if(particleB->shape[1] >= std::max(particleB->shape[0], particleB->shape[2]))
      main_axis_j = 1;
    else
      main_axis_j = 2;
    //printf("!");
  }

  const double l_i = particleA->shape[main_axis_i];
  const double l_j = particleB->shape[main_axis_j];

  const double r1_i = particleA->shape[(main_axis_i+1)%3];
  const double r2_i = particleA->shape[(main_axis_i+2)%3];
  const double r1_j = particleB->shape[(main_axis_j+1)%3];
  const double r2_j = particleB->shape[(main_axis_j+2)%3];

  double r_i = std::max(r1_i,r2_i);
  double r_j = std::max(r1_j,r2_j);
  
  if(!particleA->isCylinder)
    r_i = sqrt(r1_i*r1_i + r2_i*r2_i);
  if(!particleB->isCylinder)
    r_j = sqrt(r1_j*r1_j + r2_j*r2_j);
      
  const double radsum = r_i + r_j;

  double Pi[3], Qi[3], Pj[3], Qj[3];
  for(int k = 0; k < 3; k++) {
    Pi[k] = -particleA->rotation_matrix[k*3+main_axis_i]*l_i + particleA->center[k];
    Qi[k] =  particleA->rotation_matrix[k*3+main_axis_i]*l_i + particleA->center[k];
    Pj[k] = -particleB->rotation_matrix[k*3+main_axis_j]*l_j + particleB->center[k];
    Qj[k] =  particleB->rotation_matrix[k*3+main_axis_j]*l_j + particleB->center[k];
  }

  double C1[3], C2[3];
  MathExtraLiggghtsNonspherical::line_segments_distance(Pi, Qi, Pj, Qj, C1, C2);
  double alpha = r_i / (r_i + r_j);
  for(int k = 0; k < 3; k++)
    capsule_contact_point[k] = alpha*particleB->center[k] + (1.0-alpha)*particleA->center[k];

  return LAMMPS_NS::pointDistance(C1, C2) < radsum;
}

//calculate the contact point if no information about from the previous step available
bool calc_contact_point_if_no_previous_point_avaialable(SurfacesIntersectData & sidata, Superquadric *particleA,
    Superquadric *particleB, double *contact_point, double &fi, double &fj, LAMMPS_NS::Error *error)
{
  bool fail_flag = false;
  const double ai = particleA->shape[0];
  const double bi = particleA->shape[1];
  const double ci = particleA->shape[2];

  const double aj = particleB->shape[0];
  const double bj = particleB->shape[1];
  const double cj = particleB->shape[2];

  const double ri = cbrt(ai*bi*ci); // "effective" radius of the particle i
  const double rj = cbrt(aj*bj*cj); // "effective" radius of the particle j
  double ratio = ri / (ri + rj);
  const double ni0 = particleA->blockiness[0];
  const double ni1 = particleA->blockiness[1];
  const double nj0 = particleB->blockiness[0];
  const double nj1 = particleB->blockiness[1];

  double pre_estimation[3];
  double nmax = std::max(std::max(std::max(ni0, ni1), nj0), nj1);
  double step = 1.0;
  int N1 = (nmax - 2.0)/step + 1;
  int N = std::max(N1, 10);
  double N_inv = 1.0/static_cast<double>(N);
  for(int i = 0; i <= N; i++) { //iterate from spheres to actual shaped particles
    particleA->set_blockiness(2.0 + (ni0 - 2.0)*static_cast<double>(i) * N_inv, 2.0 + (ni1 - 2.0)*static_cast<double>(i) * N_inv);
    particleB->set_blockiness(2.0 + (nj0 - 2.0)*static_cast<double>(i) * N_inv, 2.0 + (nj1 - 2.0)*static_cast<double>(i) * N_inv);
    particleA->set_shape(ri + (ai - ri)*static_cast<double>(i) * N_inv, ri + (bi - ri)*static_cast<double>(i) * N_inv, ri + (ci - ri)*static_cast<double>(i) * N_inv);
    particleB->set_shape(rj + (aj - rj)*static_cast<double>(i) * N_inv, rj + (bj - rj)*static_cast<double>(i) * N_inv, rj + (cj - rj)*static_cast<double>(i) * N_inv);

    if(i == 0) {
      for(int k = 0; k < 3; k++)
        contact_point[k] = ratio*particleB->center[k] + (1.0 - ratio)*particleA->center[k]; //contact point solution for spheres
    }
    LAMMPS_NS::vectorCopy3D(contact_point, pre_estimation);
    MathExtraLiggghtsNonspherical::calc_contact_point(particleA, particleB,
        ratio, pre_estimation, contact_point, fi, fj, &fail_flag, error);

    particleA->set_shape(ai, bi, ci); //set actual shape parameters back
    particleB->set_shape(aj, bj, cj);
    particleA->set_blockiness(ni0, ni1);
    particleB->set_blockiness(nj0, nj1);
  }
  return fail_flag;
}

//calculate the contact point using information from previous step
bool calc_contact_point_using_prev_step(SurfacesIntersectData & sidata, Superquadric *particleA, Superquadric *particleB,
    double ratio, double dt, double *prev_step_point, double *contact_point, double &fi, double &fj, LAMMPS_NS::Error *error)
{
  bool fail_flag = false;

  #ifdef SUPERQUADRIC_ACTIVE_FLAG

  double estimation[3];
  double xci[3], xcj[3], v_rot_i[3], v_rot_j[3], vi[3], vj[3];
  LAMMPS_NS::vectorSubtract3D(prev_step_point, particleA->center, xci);
  LAMMPS_NS::vectorSubtract3D(prev_step_point, particleB->center, xcj);
  LAMMPS_NS::vectorCross3D(sidata.omega_i, xci, v_rot_i);
  LAMMPS_NS::vectorCross3D(sidata.omega_j, xcj, v_rot_j);

  for(int k = 0; k < 3; k++) {
    vi[k] = sidata.v_i[k] + v_rot_i[k]; //velocity of the contact point with respect to particle i
    vj[k] = sidata.v_j[k] + v_rot_j[k]; //velocity of the contact point with respect to particle j
  }
  double v_eff[3];
  for(int k = 0; k < 3; k++)
    v_eff[k] = (vi[k]*sidata.mi + vj[k]*sidata.mj) / (sidata.mi+sidata.mj); //average velocity of the contact point

  for(int k = 0; k < 3; k++)
    estimation[k] = prev_step_point[k] + v_eff[k]*dt; //contact point estimation

  MathExtraLiggghtsNonspherical::calc_contact_point(particleA, particleB, ratio, estimation, contact_point, fi, fj, &fail_flag, error);
  #endif

  return fail_flag;
}

//basic estimation of the overlap magnitude
void basic_overlap_algorithm(SurfacesIntersectData & sidata, Superquadric *particleA, Superquadric *particleB,
    double &alphaA, double &alphaB, const double *contact_point, double *contact_pointA, double *contact_pointB)
{
  particleA->map_point(contact_point, contact_pointA);
  particleB->map_point(contact_point, contact_pointB);
  double deltaA[3], deltaB[3];
  LAMMPS_NS::vectorSubtract3D(contact_pointA, contact_point, deltaA);
  LAMMPS_NS::vectorSubtract3D(contact_pointB, contact_point, deltaB);
  alphaA = LAMMPS_NS::vectorDot3D(deltaA, sidata.en);
  alphaB = LAMMPS_NS::vectorDot3D(deltaB, sidata.en);

  sidata.deltan = fabs(alphaA) + fabs(alphaB);

  LAMMPS_NS::vectorScalarMult3D(sidata.en, sidata.deltan, sidata.delta);
}

//more sophisticated overlap algorithm, uses the same overlap direction
double extended_overlap_algorithm(Superquadric *particleA, Superquadric *particleB, double *en,
    double *const alpha_A, double *const alpha_B,
    const double *contact_point, double *contact_pointA, double *contact_pointB, double *delta)
{
  #ifdef SUPERQUADRIC_ACTIVE_FLAG
  bool use_alpha_ellipsoid = false;
  LAMMPS_NS::vectorNegate3D(en);
  *alpha_A = particleA->surface_line_intersection(use_alpha_ellipsoid, contact_point, en, *alpha_A, contact_pointA);
  LAMMPS_NS::vectorNegate3D(en);
  *alpha_B = particleB->surface_line_intersection(use_alpha_ellipsoid, contact_point, en, *alpha_B, contact_pointB);

  LAMMPS_NS::vectorSubtract3D(contact_pointB, contact_pointA, delta);
  return fabs(*alpha_A) + fabs(*alpha_B); //overlap magnitude
  #else
  return 0.0;
  #endif
}

//common normal algorithm
double common_normal(SurfacesIntersectData & sidata, Superquadric *particleA, Superquadric *particleB, bool particles_were_in_contact_flag,
    double *const contact_point_A_history, double *const contact_point_B_history, //from previous timestep
    double *contact_point_A, double *contact_point_B) //from line-surface intersection algorithm
{
  bool fail_flag;
  double contact_point_A_temp1[3], contact_point_B_temp1[3];
  double contact_point_A_temp2[3], contact_point_B_temp2[3];

  LAMMPS_NS::vectorCopy3D(contact_point_A, contact_point_A_temp2);
  LAMMPS_NS::vectorCopy3D(contact_point_B, contact_point_B_temp2);

  if(!particles_were_in_contact_flag) {
    LAMMPS_NS::vectorCopy3D(contact_point_A_temp2, contact_point_A_temp1); //use solution from overlap model 1 as initial estimation
    LAMMPS_NS::vectorCopy3D(contact_point_B_temp2, contact_point_B_temp1);
  } else  {
    LAMMPS_NS::vectorCopy3D(contact_point_A_history, contact_point_A_temp1); //use pair of contact points from previous step as initial estimation
    LAMMPS_NS::vectorCopy3D(contact_point_B_history, contact_point_B_temp1);
  }

  double result_pointA1[3],result_pointB1[3];
  double result_pointA2[3],result_pointB2[3];

  double res1 = MathExtraLiggghtsNonspherical::minimal_distance(particleA, particleB, contact_point_A_temp1, contact_point_B_temp1,
                                                                                      result_pointA1,       result_pointB1, &fail_flag);
  double res = res1;

  if(res > 1e-12) {
    double res2 = MathExtraLiggghtsNonspherical::minimal_distance(particleA, particleB, contact_point_A_temp2, contact_point_B_temp2,
                                                                                          result_pointA2,       result_pointB2, &fail_flag);

    if(res2 > res1) {
      LAMMPS_NS::vectorCopy3D(result_pointA1, contact_point_A);
      LAMMPS_NS::vectorCopy3D(result_pointB1, contact_point_B);
      res = res1;
    } else {
      LAMMPS_NS::vectorCopy3D(result_pointA2, contact_point_A);
      LAMMPS_NS::vectorCopy3D(result_pointB2, contact_point_B);
      res = res2;
    }
  } else {
    LAMMPS_NS::vectorCopy3D(result_pointA1, contact_point_A);
    LAMMPS_NS::vectorCopy3D(result_pointB1, contact_point_B);
    res = res1;
}

  LAMMPS_NS::vectorSubtract3D(contact_point_B, contact_point_A, sidata.delta);
  LAMMPS_NS::vectorCopy3D(sidata.delta, sidata.en);  // overlap direction
  LAMMPS_NS::vectorNormalize3D(sidata.en);

  LAMMPS_NS::vectorCopy3D(contact_point_A, contact_point_A_history); //store pair of contact points in local coordinates for each particles
  LAMMPS_NS::vectorCopy3D(contact_point_B, contact_point_B_history);

  sidata.contact_point[0] = 0.5*(contact_point_A[0] + contact_point_B[0]);
  sidata.contact_point[1] = 0.5*(contact_point_A[1] + contact_point_B[1]);
  sidata.contact_point[2] = 0.5*(contact_point_A[2] + contact_point_B[2]);

  return LAMMPS_NS::vectorMag3D(sidata.delta); //overlap magnitude
  }

#ifdef LIGGGHTS_DEBUG
void printf_debug_data(Superquadric *particleA, Superquadric *particleB, double *initial_guess ,LAMMPS_NS::Error *error) //TODO: implement as warnings
{
  printf("Initial guess in Newton method: %1.13f %1.13f %1.13f\n", initial_guess[0], initial_guess[1], initial_guess[2]);
  printf("Particle 1: %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f\n",
          particleA->center[0], particleA->center[1], particleA->center[2],
          particleA->shape[0], particleA->shape[1], particleA->shape[2],
          particleA->blockiness[0], particleA->blockiness[1],
          particleA->quat[0], particleA->quat[1], particleA->quat[2], particleA->quat[3]);
  printf("Particle 2: %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f, %1.13f\n",
          particleB->center[0], particleB->center[1], particleB->center[2],
          particleB->shape[0], particleB->shape[1], particleB->shape[2],
          particleB->blockiness[0], particleB->blockiness[1],
          particleB->quat[0], particleB->quat[1], particleB->quat[2], particleB->quat[3]);
}
#endif

}
#endif
