//#**************************************************************
//#
//# filename:             linalgeig.cpp
//#
//# project:              APART
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						April 2006
//# description:          Eigenvalue solver from TNT library
//#                       
//# remarks:						  
//#
//# This file is part of the program package HOTINT and underlies the stipulations
//# of the scientific or license agreement. It is therefore emphasized not to copy
//# or redistribute this file. The use of this file is only permitted for academic or scholar
//# research. It is forbidden to use any part of this code for military applications!
//# The Developer does not assume any liability for this code or for results obtained
//# within its use, nor shall he be liable for any direct, indirect, or other damage
//# resulting from use of this code. 
//#
//# bug reports are welcome!!!
//# WWW: http://tmech.mechatronik.uni-linz.ac.at/staff/gerstmayr/gerstmayr.html
//# email: jg@jku.at
//#**************************************************************
//#include "stdafx.h"
#include "ioincludes.h"
#include "windows.h" //for shell execute
//#include "shellapi.h"

#include <assert.h>
#include <memory.h>

#include <math.h>

#include "mbs_interface.h"
#include "lapack_routines.h"
#include "linalgeig.h"
#include "femathhelperfunctions.h"

void Eigenvalue::Tred2() 
{

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}

	// Householder reduction to tridiagonal form.

	for (int i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} else {

			// Generate Householder vector.

			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for (int j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for (int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.

	for (int i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for (int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (int k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}
	for (int j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.
void Eigenvalue::TQL2 () 
{

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0,-52.0);
	for (int l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = Maximum(tst1,fabs(d[l]) + fabs(e[l]));
		int m = l;

		// Original while-loop from Java code
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}


		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = _hypot(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for (int i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for (int i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = _hypot(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for (int k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.

	for (int i = 0; i < n-1; i++) {
		int k = i;
		double p = d[i];
		for (int j = i+1; j < n; j++) {
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

// Nonsymmetric reduction to Hessenberg form.

void Eigenvalue::OrtHes() 
{

	//  This is derived from the Algol procedures orthes and ortran,
	//  by Martin and Wilkinson, Handbook for Auto. Comp.,
	//  Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutines in EISPACK.

	int low = 0;
	int high = n-1;

	for (int m = low+1; m <= high-1; m++) {

		// Scale column.

		double scale = 0.0;
		for (int i = m; i <= high; i++) {
			scale = scale + fabs(H[i][m-1]);
		}
		if (scale != 0.0) {

			// Compute Householder transformation.

			double h = 0.0;
			for (int i = high; i >= m; i--) {
				ort[i] = H[i][m-1]/scale;
				h += ort[i] * ort[i];
			}
			double g = sqrt(h);
			if (ort[m] > 0) {
				g = -g;
			}
			h = h - ort[m] * g;
			ort[m] = ort[m] - g;

			// Apply Householder similarity transformation
			// H = (I-u*u'/h)*H*(I-u*u')/h)

			for (int j = m; j < n; j++) {
				double f = 0.0;
				for (int i = high; i >= m; i--) {
					f += ort[i]*H[i][j];
				}
				f = f/h;
				for (int i = m; i <= high; i++) {
					H[i][j] -= f*ort[i];
				}
			}

			for (int i = 0; i <= high; i++) {
				double f = 0.0;
				for (int j = high; j >= m; j--) {
					f += ort[j]*H[i][j];
				}
				f = f/h;
				for (int j = m; j <= high; j++) {
					H[i][j] -= f*ort[j];
				}
			}
			ort[m] = scale*ort[m];
			H[m][m-1] = scale*g;
		}
	}

	// Accumulate transformations (Algol's ortran).

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			V[i][j] = (i == j ? 1.0 : 0.0);
		}
	}

	for (int m = high-1; m >= low+1; m--) {
		if (H[m][m-1] != 0.0) {
			for (int i = m+1; i <= high; i++) {
				ort[i] = H[i][m-1];
			}
			for (int j = m; j <= high; j++) {
				double g = 0.0;
				for (int i = m; i <= high; i++) {
					g += ort[i] * V[i][j];
				}
				// Double division avoids possible underflow
				g = (g / ort[m]) / H[m][m-1];
				for (int i = m; i <= high; i++) {
					V[i][j] += g * ort[i];
				}
			}
		}
	}
}


// Complex scalar division.

void Eigenvalue::Cdiv(double xr, double xi, double yr, double yi, double& cdivr, double& cdivi)
{
	double r,d;
	if (fabs(yr) > fabs(yi)) {
		r = yi/yr;
		d = yr + r*yi;
		cdivr = (xr + r*xi)/d;
		cdivi = (xi - r*xr)/d;
	} else {
		r = yr/yi;
		d = yi + r*yr;
		cdivr = (r*xr + xi)/d;
		cdivi = (r*xi - xr)/d;
	}
}


// Nonsymmetric reduction from Hessenberg to double Schur form.

void Eigenvalue::HQR2 () 
{

	//  This is derived from the Algol procedure hqr2,
	//  by Martin and Wilkinson, Handbook for Auto. Comp.,
	//  Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	// Initialize

	int nn = this->n;
	int n = nn-1;
	int low = 0;
	int high = nn-1;
	double eps = pow(2.0,-52.0);
	double exshift = 0.0;
	double p=0,q=0,r=0,s=0,z=0,t,w,x,y;
	double cdivr, cdivi;

	// Store roots isolated by balanc and compute matrix norm

	double norm = 0.0;
	for (int i = 0; i < nn; i++) {
		if ((i < low) || (i > high)) {
			d[i] = H[i][i];
			e[i] = 0.0;
		}
		for (int j = Maximum(i-1,0); j < nn; j++) {
			norm = norm + fabs(H[i][j]);
		}
	}

	// Outer loop over eigenvalue index

	int iter = 0;
	while (n >= low) {

		// Look for single small sub-diagonal element

		int l = n;
		while (l > low) {
			s = fabs(H[l-1][l-1]) + fabs(H[l][l]);
			if (s == 0.0) {
				s = norm;
			}
			if (fabs(H[l][l-1]) < eps * s) {
				break;
			}
			l--;
		}

		// Check for convergence
		// One root found

		if (l == n) {
			H[n][n] = H[n][n] + exshift;
			d[n] = H[n][n];
			e[n] = 0.0;
			n--;
			iter = 0;

			// Two roots found

		} else if (l == n-1) {
			w = H[n][n-1] * H[n-1][n];
			p = (H[n-1][n-1] - H[n][n]) / 2.0;
			q = p * p + w;
			z = sqrt(fabs(q));
			H[n][n] = H[n][n] + exshift;
			H[n-1][n-1] = H[n-1][n-1] + exshift;
			x = H[n][n];

			// double pair

			if (q >= 0) {
				if (p >= 0) {
					z = p + z;
				} else {
					z = p - z;
				}
				d[n-1] = x + z;
				d[n] = d[n-1];
				if (z != 0.0) {
					d[n] = x - w / z;
				}
				e[n-1] = 0.0;
				e[n] = 0.0;
				x = H[n][n-1];
				s = fabs(x) + fabs(z);
				p = x / s;
				q = z / s;
				r = sqrt(p * p+q * q);
				p = p / r;
				q = q / r;

				// Row modification

				for (int j = n-1; j < nn; j++) {
					z = H[n-1][j];
					H[n-1][j] = q * z + p * H[n][j];
					H[n][j] = q * H[n][j] - p * z;
				}

				// Column modification

				for (int i = 0; i <= n; i++) {
					z = H[i][n-1];
					H[i][n-1] = q * z + p * H[i][n];
					H[i][n] = q * H[i][n] - p * z;
				}

				// Accumulate transformations

				for (int i = low; i <= high; i++) {
					z = V[i][n-1];
					V[i][n-1] = q * z + p * V[i][n];
					V[i][n] = q * V[i][n] - p * z;
				}

				// Complex pair

			} else {
				d[n-1] = x + p;
				d[n] = x + p;
				e[n-1] = z;
				e[n] = -z;
			}
			n = n - 2;
			iter = 0;

			// No convergence yet

		} else {

			// Form shift

			x = H[n][n];
			y = 0.0;
			w = 0.0;
			if (l < n) {
				y = H[n-1][n-1];
				w = H[n][n-1] * H[n-1][n];
			}

			// Wilkinson's original ad hoc shift

			if (iter == 10) {
				exshift += x;
				for (int i = low; i <= n; i++) {
					H[i][i] -= x;
				}
				s = fabs(H[n][n-1]) + fabs(H[n-1][n-2]);
				x = y = 0.75 * s;
				w = -0.4375 * s * s;
			}

			// MATLAB's new ad hoc shift

			if (iter == 30) {
				s = (y - x) / 2.0;
				s = s * s + w;
				if (s > 0) {
					s = sqrt(s);
					if (y < x) {
						s = -s;
					}
					s = x - w / ((y - x) / 2.0 + s);
					for (int i = low; i <= n; i++) {
						H[i][i] -= s;
					}
					exshift += s;
					x = y = w = 0.964;
				}
			}

			iter = iter + 1;   // (Could check iteration count here.)

			// Look for two consecutive small sub-diagonal elements

			int m = n-2;
			while (m >= l) {
				z = H[m][m];
				r = x - z;
				s = y - z;
				p = (r * s - w) / H[m+1][m] + H[m][m+1];
				q = H[m+1][m+1] - z - r - s;
				r = H[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p = p / s;
				q = q / s;
				r = r / s;
				if (m == l) {
					break;
				}
				if (fabs(H[m][m-1]) * (fabs(q) + fabs(r)) <
					eps * (fabs(p) * (fabs(H[m-1][m-1]) + fabs(z) +
					fabs(H[m+1][m+1])))) {
						break;
				}
				m--;
			}

			for (int i = m+2; i <= n; i++) {
				H[i][i-2] = 0.0;
				if (i > m+2) {
					H[i][i-3] = 0.0;
				}
			}

			// Double QR step involving rows l:n and columns m:n

			for (int k = m; k <= n-1; k++) {
				int notlast = (k != n-1);
				if (k != m) {
					p = H[k][k-1];
					q = H[k+1][k-1];
					r = (notlast ? H[k+2][k-1] : 0.0);
					x = fabs(p) + fabs(q) + fabs(r);
					if (x != 0.0) {
						p = p / x;
						q = q / x;
						r = r / x;
					}
				}
				if (x == 0.0) {
					break;
				}
				s = sqrt(p * p + q * q + r * r);
				if (p < 0) {
					s = -s;
				}
				if (s != 0) {
					if (k != m) {
						H[k][k-1] = -s * x;
					} else if (l != m) {
						H[k][k-1] = -H[k][k-1];
					}
					p = p + s;
					x = p / s;
					y = q / s;
					z = r / s;
					q = q / p;
					r = r / p;

					// Row modification

					for (int j = k; j < nn; j++) {
						p = H[k][j] + q * H[k+1][j];
						if (notlast) {
							p = p + r * H[k+2][j];
							H[k+2][j] = H[k+2][j] - p * z;
						}
						H[k][j] = H[k][j] - p * x;
						H[k+1][j] = H[k+1][j] - p * y;
					}

					// Column modification

					for (int i = 0; i <= Minimum(n,k+3); i++) {
						p = x * H[i][k] + y * H[i][k+1];
						if (notlast) {
							p = p + z * H[i][k+2];
							H[i][k+2] = H[i][k+2] - p * r;
						}
						H[i][k] = H[i][k] - p;
						H[i][k+1] = H[i][k+1] - p * q;
					}

					// Accumulate transformations

					for (int i = low; i <= high; i++) {
						p = x * V[i][k] + y * V[i][k+1];
						if (notlast) {
							p = p + z * V[i][k+2];
							V[i][k+2] = V[i][k+2] - p * r;
						}
						V[i][k] = V[i][k] - p;
						V[i][k+1] = V[i][k+1] - p * q;
					}
				}  // (s != 0)
			}  // k loop
		}  // check convergence
	}  // while (n >= low)

	// Backsubstitute to find vectors of upper triangular form

	if (norm == 0.0) {
		return;
	}

	for (n = nn-1; n >= 0; n--) {
		p = d[n];
		q = e[n];

		// double vector

		if (q == 0) {
			int l = n;
			H[n][n] = 1.0;
			for (int i = n-1; i >= 0; i--) {
				w = H[i][i] - p;
				r = 0.0;
				for (int j = l; j <= n; j++) {
					r = r + H[i][j] * H[j][n];
				}
				if (e[i] < 0.0) {
					z = w;
					s = r;
				} else {
					l = i;
					if (e[i] == 0.0) {
						if (w != 0.0) {
							H[i][n] = -r / w;
						} else {
							H[i][n] = -r / (eps * norm);
						}

						// Solve double equations

					} else {
						x = H[i][i+1];
						y = H[i+1][i];
						q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
						t = (x * s - z * r) / q;
						H[i][n] = t;
						if (fabs(x) > fabs(z)) {
							H[i+1][n] = (-r - w * t) / x;
						} else {
							H[i+1][n] = (-s - y * t) / z;
						}
					}

					// Overflow control

					t = fabs(H[i][n]);
					if ((eps * t) * t > 1) {
						for (int j = i; j <= n; j++) {
							H[j][n] = H[j][n] / t;
						}
					}
				}
			}

			// Complex vector

		} else if (q < 0) {
			int l = n-1;

			// Last vector component imaginary so matrix is triangular

			if (fabs(H[n][n-1]) > fabs(H[n-1][n])) {
				H[n-1][n-1] = q / H[n][n-1];
				H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
			} else {
				Cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q, cdivr, cdivi);
				H[n-1][n-1] = cdivr;
				H[n-1][n] = cdivi;
			}
			H[n][n-1] = 0.0;
			H[n][n] = 1.0;
			for (int i = n-2; i >= 0; i--) {
				double ra,sa,vr,vi;
				ra = 0.0;
				sa = 0.0;
				for (int j = l; j <= n; j++) {
					ra = ra + H[i][j] * H[j][n-1];
					sa = sa + H[i][j] * H[j][n];
				}
				w = H[i][i] - p;

				if (e[i] < 0.0) {
					z = w;
					r = ra;
					s = sa;
				} else {
					l = i;
					if (e[i] == 0) {
						Cdiv(-ra,-sa,w,q, cdivr, cdivi);
						H[i][n-1] = cdivr;
						H[i][n] = cdivi;
					} else {

						// Solve complex equations

						x = H[i][i+1];
						y = H[i+1][i];
						vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
						vi = (d[i] - p) * 2.0 * q;
						if ((vr == 0.0) && (vi == 0.0)) {
							vr = eps * norm * (fabs(w) + fabs(q) +
								fabs(x) + fabs(y) + fabs(z));
						}
						Cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi, cdivr, cdivi);
						H[i][n-1] = cdivr;
						H[i][n] = cdivi;
						if (fabs(x) > (fabs(z) + fabs(q))) {
							H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
							H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
						} else {
							Cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q, cdivr, cdivi);
							H[i+1][n-1] = cdivr;
							H[i+1][n] = cdivi;
						}
					}

					// Overflow control

					t = Maximum(fabs(H[i][n-1]),fabs(H[i][n]));
					if ((eps * t) * t > 1) {
						for (int j = i; j <= n; j++) {
							H[j][n-1] = H[j][n-1] / t;
							H[j][n] = H[j][n] / t;
						}
					}
				}
			}
		}
	}

	// Vectors of isolated roots

	for (int i = 0; i < nn; i++) {
		if (i < low || i > high) {
			for (int j = i; j < nn; j++) {
				V[i][j] = H[i][j];
			}
		}
	}

	// Back transformation to get eigenvectors of original matrix

	for (int j = nn-1; j >= low; j--) {
		for (int i = low; i <= high; i++) {
			z = 0.0;
			for (int k = low; k <= Minimum(j,high); k++) {
				z = z + V[i][k] * H[k][j];
			}
			V[i][j] = z;
		}
	}
}


/** Check for symmetry, then construct the eigenvalue decomposition
@param A    Square double (non-complex) matrix
*/
Eigenvalue::Eigenvalue(const Matrix &A) 
{
#ifndef __QUICKMATH
	assert(A.IsSquare());
#endif
	n = A.Getcols();
	V = Matrix(n,n);
	d = Vector(n);
	e = Vector(n);

	issymmetric = 1;
	for (int j = 0; (j < n) && issymmetric; j++) {
		for (int i = 0; (i < n) && issymmetric; i++) {
			issymmetric = (A[i][j] == A[j][i]);
		}
	}

	if (issymmetric) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V[i][j] = A[i][j];
			}
		}

		// Tridiagonalize.
		Tred2();

		// Diagonalize.
		TQL2();

	} else {
		H = Matrix(n,n);
		ort = Vector(n);

		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				H[i][j] = A[i][j];
			}
		}

		// Reduce to Hessenberg form.
		OrtHes();

		// Reduce Hessenberg to double Schur form.
		HQR2();
	}
}



/** 
Computes the block diagonal eigenvalue matrix.
If the original matrix A is not symmetric, then the eigenvalue 
matrix D is block diagonal with the real eigenvalues in 1-by-1 
blocks and any complex eigenvalues,
a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
eigenvalues look like
<pre>

u + iv     .        .          .      .    .
.      u - iv     .          .      .    .
.        .      a + ib       .      .    .
.        .        .        a - ib   .    .
.        .        .          .      x    .
.        .        .          .      .    y
</pre>
then D looks like
<pre>

u        v        .          .      .    .
-v        u        .          .      .    . 
.        .        a          b      .    .
.        .       -b          a      .    .
.        .        .          .      x    .
.        .        .          .      .    y
</pre>
This keeps V a double matrix in both symmetric and non-symmetric
cases, and A*V = V*D.

@param D: upon return, the matrix is filled with the block diagonal 
eigenvalue matrix.

*/

void Eigenvalue::GetD(Matrix &D) 
{
	D = Matrix(n,n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			D[i][j] = 0.0;
		}
		D[i][i] = d[i];
		if (e[i] > 0) {
			D[i][i+1] = e[i];
		} else if (e[i] < 0) {
			D[i][i-1] = e[i];
		}
	}
}

// compute nev smallest eigenvalues of the generalized eigenvalue system 
// A u = lambda M u
// by LOBPCG method (Knyazev)
// Vector eigenvalues of length nev
//   and	
// Array<Vector*> of length nev, vectors of matrix dimensions
// have to be provided!!
// returns 1 if computed correctly, else 0
int SparseEigenvalueSolver::ComputeEigenModes(int nev, Vector& eigenvalues, Matrix& eigmodes, TArray<int>& unconstraineddofs, int nzeromodes)
{
	int symmetric_lapack_evsolver = 1;
	double time = GetClockTime();
	int writetofile = 1;
	// problem size
	// matrices A, M need not be provided (NULL), thus take length of eigenmode vector as size
	int size = unconstraineddofs.Length();
	int N = 3*nev; // useful, number of x, w, p
	if ( N > size )
	{
		mbs->UO() << "Maximum number of eigenmodes which can be computed with LOBPCG method is Matrix size/3 = " << size/3 << "!\n";
		return 0;
	}

	srand(10); //initialize always with same number->always same result!!!

	ofstream outfile("..\\..\\output\\output_LOBPCG.txt");
	if (writetofile)
	{
		outfile << "LOBPCG solver, finding " << nev << " smallest eigenvalues\n";
		outfile << "Eigenvalue problem size " << size << "\n";
	}

	// Matrices where constrained degrees of freedom are removed
	// if no constrained dofs, then Kmat and Mmat are used
	// if unconstrained dofs, new matrices are allocated, which have to be deleted!
	SparseMatrix *Ared, *Mred;
	int allocateA = 0;
	int allocateM = 0;
	if(!A)
	{
		Ared = A;
		Mred = M;
	}	
	else if (size != A->Getrows())
	{
		Ared = new SparseMatrix();
		Ared->CopyFrom(*A,1,1,size,size, unconstraineddofs);
		allocateA = 1;
		if (M)
		{
			Mred = new SparseMatrix();
			Mred->CopyFrom(*M,1,1,size,size, unconstraineddofs);
			allocateM = 1;
		}
		else
		{
			Mred = 0;
		}
	}
	else
	{
		Ared = A;
		Mred = M;
	}

	// store eigenvalues from last iterative step
	// corresponds to vector mu in Knyazev
	Vector eigenvals_old(nev);
	eigenvals_old.SetAll(0.);

	// Vectors p_i, w_i, needed in iteration 
	TArray<Vector*> p, w, pold, eigenmodes, eigenmodes_old;
	p.SetLen(nev);
	w.SetLen(nev);
	pold.SetLen(nev);
	eigenmodes_old.SetLen(nev);
	for (int m=1; m<=nev; m++)
	{
		p(m) = new Vector(size);
		w(m) = new Vector(size);
		pold(m) = new Vector(size);
		eigenmodes_old(m) = new Vector(size);
		eigenmodes(m) = new Vector(size);
	}

	// help vector A*u, M*u
	Vector Au(size);
	Vector Mu(size);
	// eigenvalues: in nonsymmetric case, eigenvalues are Complex(ev_small, ev_smalli) / ev_beta
	// TArray is used for sorting reasons (Quicksort)
	TArray<double> ev_small(N);
	TArray<double> ev_smalli(N);
	TArray<double> ev_beta(N);
	// resorted eigenvalue indices
	TArray<int> sortref(N);
	for (int i=1; i<=N; i++)
		sortref(i) = i;
	// matrices for small projected eigenvalue problems
	Matrix Asmall(N,N);
	Matrix EVmat(N,N);
	Matrix Msmall(N,N);
	Vector work(8*N+16);

	// Initialize vectors:
	// eigenmodes are initialized randomly if flag is set
	if (InitializeRandom())
	{
		for (int i=1; i<=nev; i++)
		{
			eigenmodes(i)->SetRandom();

		}
	}
	// p, w are initialized with zeros
	for (int i=1; i<=nev; i++)
	{
		p(i)->FillWithZeros();
		w(i)->FillWithZeros();
	}



	// first iteration: p = 0, N is number of x, w
	// afterwards, N is reset to 3*nev = number of x,w,p
	N = 2*nev;

	//Iteration
	int step = 1;
	double error = 1;
	double error_zeromodes = 1;
	double error_init = 0.1*error/precision;

	SparseInverse *invA = 0;
	if (use_preconditioner)
	{
		if (Mred)
			Ared->AddSubmatrix(*Mred,1,1,1,1,size,size,lambda_precond);
		else
		{
			for (int i=1; i<=size; i++)
				(*Ared)(i,i) += lambda_precond;
		}
		invA = new SparseInverse(*Ared);
		invA->Factorize();
		if (Mred)
			Ared->AddSubmatrix(*Mred,1,1,1,1,size,size,-lambda_precond);
		else
		{
			for (int i=1; i<=size; i++)
				(*Ared)(i,i) -= lambda_precond;
		}
	}

	int LAPACKINFO = 0;
	bool stopcalculation = mbs->StopCalculation() && CanBeStopped();

	if (mbs->UO().GetGlobalMessageLevel() >= UO_LVL_ext)
		mbs->UO() << "starting iteration\n";
	
	int convergedcorrectly=0;
	while(step <= MaxSteps() && !convergedcorrectly && !LAPACKINFO && !stopcalculation)
	{

		error = 0;
		error_zeromodes = 0;
		for(int m=1; m<=nev; m++)
		{
			// compute A*x(m) in w(m), M*x(m) in Mu
			ApplyMat(Mred,*eigenmodes(m), Mu);
			ApplyMat(Ared,*eigenmodes(m), *w(m));

			// -- Step 3. --  In first iterative step, compute eigenvalue approximations mu for initial eigenmode guess
			double xTAx = *eigenmodes(m) * (*w(m));
			double xTMx = *eigenmodes(m) * Mu;
			eigenvalues(m) = xTAx/xTMx; 

			// -- Step 4. --  Compute residual vector in w(m)
			// r = Au - eigenvalues(m) * Mu;
			w(m)->MultAdd(-eigenvalues(m), Mu);

			// -- Step 5. -- Preconditioning	
			// w = Tr
			if (use_preconditioner)
				invA->Apply(*w(m));

			if (m <= nzeromodes)
			{
				error_zeromodes += fabs(eigenvalues(m));
			}
			else
			{
				error += w(m)->GetNorm();
			}
			// normalize residual direction
			*w(m) *= 1./w(m)->GetNorm();
		}

		if (step==1)
		{
			if (error)
				error_init = error;
			else
				error_init = 1;
		}
		//if (step%1==0 && (mbs->UO().GetGlobalMessageLevel() >= UO_LVL_ext)) // old code --> not very clever
		int dont_show_message = step % 10;
		if ((!dont_show_message) || (mbs->UO().GetGlobalMessageLevel() >= UO_LVL_ext))	//$ DR 2013-03-29
			mbs->UO(UO_LVL_sim) << "Step " << step << ": relative error " << error/error_init << ", error_zeromodes = " << error_zeromodes << "\n";

		// -- Step 6. --  Eigenvalue computation for small matrix on test space [x1..xn, w1..wn, p1..pn]
		// -- -- build small block matrices
		Asmall.SetSize(N,N);
		Msmall.SetSize(N,N);
		Asmall.SetAll(0.);
		Msmall.SetAll(0.);
		for (int i=1; i<=nev; i++)
		{
			// XX-block
			ApplyMat(Ared,*eigenmodes(i), Au);
			ApplyMat(Mred,*eigenmodes(i), Mu);
			for (int k=i; k<=nev; k++)
			{
				Asmall(i,k) += *eigenmodes(k)*Au;
				Msmall(i,k) += *eigenmodes(k)*Mu;
			}
			// XW-block
			for (int k=1; k<=nev; k++)
			{
				Asmall(i,nev+k) += *w(k)*Au;
				Msmall(i,nev+k) += *w(k)*Mu;
			}
			// XP-block
			if (N==3*nev)
			{
				for (int k=1; k<=nev; k++)
				{
					Asmall(i,2*nev+k) += *p(k)*Au;
					Msmall(i,2*nev+k) += *p(k)*Mu;
				}
			}
			// WW-block
			ApplyMat(Ared,*w(i), Au);
			ApplyMat(Mred,*w(i), Mu);
			for (int k=i; k<=nev; k++)
			{
				Asmall(nev+i,nev+k) += *w(k)*Au;
				Msmall(nev+i,nev+k) += *w(k)*Mu;
			}
			if (N==3*nev)
			{
				// WP-block
				for (int k=1; k<=nev; k++)
				{
					Asmall(nev+i,2*nev+k) += *p(k)*Au;
					Msmall(nev+i,2*nev+k) += *p(k)*Mu;
				}
				// PP-block
				ApplyMat(Ared,*p(i), Au);
				ApplyMat(Mred,*p(i), Mu);
				for (int k=i; k<=nev; k++)
				{
					Asmall(2*nev+i,2*nev+k) += *p(k)*Au;
					Msmall(2*nev+i,2*nev+k) += *p(k)*Mu;
				}
			}
		}
		//// for LAPACK EV-solver, only upper right triangular matrix is needed
		
		if (!symmetric_lapack_evsolver)
		{
		for (int i=1; i<=N; i++)
			for (int j=1; j<i; j++)
			{
				Asmall(i,j) = Asmall(j,i);
				Msmall(i,j) = Msmall(j,i);
			}
		}

		EVmat.SetSize(N,N);
		// -- -- compute block matrix eigenvalues
		// solve gen. eigenvalue problem, eigenvectors are stored rowwise in Asmall
		// eigenvalues are sorted by size, smallest EV first
		if (symmetric_lapack_evsolver)
		{
			LAPACKINFO = LapackGenEVPSPD(N, &Asmall(1,1), &Msmall(1,1), &ev_small(1), &work(1), work.Length());
		}
		else
		{
			LAPACKINFO = LapackGenEVP(N, &Asmall(1,1), &Msmall(1,1), &ev_small(1), &ev_smalli(1), &ev_beta(1), &EVmat(1,1), &work(1), work.Length());
			// sort eigenvalues by size
			Asmall = EVmat;
			for (int i=1; i<=N; i++)
			{
				if (fabs(ev_smalli(i)) > 1e-14)
				{
					mbs->UO() << "LapackGenEVP: complex eigenvalue obtained: lambda(" << i << ") = " << ev_small(i) / ev_beta(i) <<" + "<< ev_smalli(i)/ev_beta(i) << "\n";
					mbs->UO() << "  Computation is stopped!\n";
					LAPACKINFO = -100;
				}
				if (fabs(ev_beta(i)) < 1e-14)
				{
					mbs->UO() << "LapackGenEVP: Matrix M (" << N << "x" << N << ")  of intermediate eigenvalue problem is nearly singular, computation is stopped!\n";
					LAPACKINFO = -200;
			}
				ev_beta(i) = ev_small(i) / ev_beta(i);
			}

			sortref.SetLen(N);
			for (int i=1; i<=N; i++)
				sortref(i) = i;
			QuicksortDouble(ev_beta, sortref);
		}

		if(LAPACKINFO)
		{
			//mbs->UO().InstantMessageText("LapackGenEVP: eigenvalue problem not solved correctly!");
			if (LAPACKINFO < 0 && LAPACKINFO > -99)
			{
				mbs->UO() << "Lapack: Incorrect input value at parameter " << -LAPACKINFO << ". Computation is stopped!\n";
			}
			else if (LAPACKINFO <= N && symmetric_lapack_evsolver)
			{
				mbs->UO() << "Lapack: " << LAPACKINFO << " off-diagonal elements did not converge to zero! Computation is stopped!\n";
			}
			else if (LAPACKINFO <= 2*N && symmetric_lapack_evsolver)
			{
				mbs->UO() << "Lapack: " << LAPACKINFO-N << "th leading minor of intermediate matrix M (" << N << "x" << N << ") is not positive definite! Computation is stopped!\n";
			}
			else if (LAPACKINFO <= N && !symmetric_lapack_evsolver)
			{
				mbs->UO() << "Lapack: QZ iteration failed! Computation is stopped!\n";
			}
			else if (LAPACKINFO > N && !symmetric_lapack_evsolver)
			{
				mbs->UO() << "Lapack: Internal error!\n";
			}
			mbs->UO() << "LapackGenEVP: eigenvalue problem not solved correctly ";
			mbs->UO() << "in step " << step << ", info = " << LAPACKINFO << "!\n";
			break;
		}



		//for (int i=1; i<=nev; i++)
		//	eigenvalues(i) = ev_small(N-1-i);
		// normalize eigenmodes
		double sign = 1;
		int foundsign = 0;
		for(int j=1; j <= N; j++)
		{
			double invnorm = 0;
			Vector Aj;
			Aj.LinkWith(&Asmall(j,1),N);
			invnorm = Aj.GetNorm();
			if (invnorm == 0) {invnorm = 1;}
			invnorm = 1./invnorm;
			Aj *= invnorm;
		}


		// -- Step 7. --  compute new approximation for eigenmode
		// and
		// -- Step 8. --  compute new update direction p
		for (int i=1; i<=nev; i++)
		{
			*pold(i) = *p(i);
			*eigenmodes_old(i) = *eigenmodes(i);
			p(i)->SetAll(0.);
			eigenmodes(i)->SetAll(0);
		}

		for (int j=1; j<=nev; j++)
		{
			for (int k=1; k<=nev; k++)
			{
				int J = sortref(j);
				eigenmodes(j)->MultAdd(Asmall(J,k), *eigenmodes_old(k));
				eigenmodes(j)->MultAdd(Asmall(J,nev+k), *w(k));
				p(j)->MultAdd(Asmall(J,nev+k), *w(k));
				if (N==3*nev)
				{
					eigenmodes(j)->MultAdd(Asmall(J,2*nev+k),*pold(k));
					p(j)->MultAdd(Asmall(J,2*nev+k), *pold(k));
				}
			}
			double invnorme = 1./eigenmodes(j)->GetNorm();

			// orthogonalize search directions
			for (int k=1; k<j; k++)
			{
				double pjpk = (*p(j))*(*p(k));
				p(j)->MultAdd(-pjpk, *p(k));
			}
			double invnormp = 1./p(j)->GetNorm();
			*eigenmodes(j) *= invnorme;
			*p(j) *= invnormp;
		}

		eigenvals_old = eigenvalues;
		if (writetofile && (step==1 || step%10==0) )
		{
			outfile << "Step " << step << ":\n";
			outfile << "** eigenvalues " << eigenvals_old << "\n";
			outfile << "** error " << error/error_init << "\n";
			outfile.flush();
		}

		N = 3*nev;


		step++;
		stopcalculation = mbs->StopCalculation() && CanBeStopped(); 

		if (error/error_init < precision && error_zeromodes < precision)
			convergedcorrectly = 1;
	}// end Iteration

	time = GetClockTime() - time;
	mbs->UO() << step << " iterations, time taken " << time << "s\n";

	for(int j=1; j <= nev; j++)//number of eigenvectors
	{
		for(int i=1; i <= size; i++)//length of eigenvectors
		{	
			eigmodes(j, unconstraineddofs(i)) = (*eigenmodes(j))(i);
		}
	}

	// delete p, w;
	for (int m=1; m<=nev; m++)
	{
		delete p(m);
		delete w(m);
		delete pold(m);
		delete eigenmodes_old(m);
		delete eigenmodes(m);
	}

	if (allocateA)
		delete Ared;
	if (allocateM)
		delete Mred;
	delete invA;

	return convergedcorrectly;
	
	// send message to close eigenvalue computation window
	

}



// compute generalized eigenmodes in Matlab
// neig smallest eigenvalues are computed iteratively
// eigenmodes are stored ROWwise in Matrix eigmodes
// eigval and eigmodes have to be provided with appropriate size!
//returns 1 if computed correctly, 0 else
int SparseEigenvalueSolver::ComputeEigenModesMatlab(const mystr& pathMatlab, int neig, 
														Vector& eigval, Matrix& eigmodes, TArray<int>& unconstraineddofs,
														int nzeromodes)
{
	// Matrices where constrained degrees of freedom are removed
	// if no constrained dofs, then Kmat and Mmat are used
	// if unconstrained dofs, new matrices are allocated, which have to be deleted!
	SparseMatrix *Kred, *Mred;
	int allocateKM = 0;
	int size = unconstraineddofs.Length();
	if (size != A->Getrows())
	{
		Kred = new SparseMatrix();
		Mred = new SparseMatrix();
		Kred->CopyFrom(*A,1,1,size,size, unconstraineddofs);
		Mred->CopyFrom(*M,1,1,size,size, unconstraineddofs);
		allocateKM = 1;
	}
	else
	{
		Kred = A;
		Mred = M;
	}

	ofstream Koutput(pathMatlab+mystr("/Kmat.dat"));//generate files in path from above
	ofstream Moutput(pathMatlab+mystr("/Mmat.dat"));
	ofstream wait(pathMatlab+mystr("/wait.txt"));
	wait << "0\n";//0 in the wait-File means that the derivation is not finished
	ofstream matout(pathMatlab+mystr("/MatlabEigenvalueSolve.m"));

	// print matrices
	Kred->PrintToMatlabSparse(Koutput);
	Koutput << flush;

	Mred->PrintToMatlabSparse(Moutput);
	Moutput << flush;

	// write file which matlab will run
	matout << "load Mmat.dat\n";
	matout << "Mmat_sparse=spconvert(Mmat);\n";//saves Mmat in a Matlab-Sparse-Format
	matout << "\n";

	matout << "load Kmat.dat\n";
	matout << "Kmat_sparse=spconvert(Kmat);\n";
	matout << "\n";

	matout << "neig=" << neig << "\n";
	matout << "maxit=" << MaxSteps() << "\n";
	matout << "tol=" << precision << "\n";
	matout << "numberzeromodes=" << nzeromodes << "\n";
	mbs->UO() << "maxit=" << MaxSteps() << "\n";
	mbs->UO() << "tol=" << precision << "\n";
	mbs->UO() << "numberzeromodes=" << nzeromodes << "\n";

	matout << "compute_eigensystem(Kmat_sparse,Mmat_sparse,neig,maxit,tol,numberzeromodes);\n";//calculates Eigensystem and writes values in 
	//matout << "compute_eigensystem(Kmat_sparse,Mmat_sparse,ConstrDOF,neig,maxit,tol,numberzeromodes);\n";//calculates Eigensystem and writes values in "Eigensystem.txt"
	//matout << "exit\n"; //closes Matlab after computation automatically
	matout << flush;


	mbs->UO() << "compute eigenvectors in matlab ... \n";
	ShellExecute(NULL, "open", "matlab"," /r MatlabEigenvalueSolve", pathMatlab.c_str(), SW_SHOW);

	//1 in wait-File -> computation is finished
	int goon = 0;
	bool stopcalculation = mbs->StopCalculation() && CanBeStopped();

	while (!goon && !stopcalculation)
	{
		{
			Sleep(1000);
			ifstream wait(pathMatlab+mystr("/wait.txt"));
			int test;
			wait >> test;
			if (test == 1) goon = 1;
			stopcalculation = mbs->StopCalculation() && CanBeStopped();
			//UO() << "not ready\n";
		}
	}
	mbs->UO() << "Eigenvectors computed in MATLAB!!!\n";



	//eigenfrequencies and eigenvectors are written to "Eigensystem.txt":
	ifstream matin(pathMatlab+mystr("/Eigensystem.txt"));
	mbs->UO() << "read eigenvectors from matlab ...\n";

	int conv;
	// read eigenmodes and values from file, if computation in matlab was finished
	if (stopcalculation)
	{
		conv = 1;
		eigval.SetAll(0);
		eigmodes.SetAll(0);
	}
	else
	{
		for (int i=1; i <= neig; i++)
		{
			matin >> eigval(i);//eigenfrequencies
		}


		for(int i=1; i <= size; i++)//length of eigenvectors
		{	
			for(int j=1; j <= neig; j++)//number of eigenvectors
			{
				matin >> eigmodes(j, unconstraineddofs(i));
				//EigenvektorMatrix(i,j)*=0.4; //deformation scaling factor = 0.4
			}
		}

		//Message for user, if convergence is not reached:
		ifstream ConvergenceInfo(pathMatlab+mystr("/Convergence.txt"));
		ConvergenceInfo >> conv;
	}
	//UO() << "convflag="<<conv<<"\n";
	if (allocateKM)
	{
		delete Kred;
		delete Mred;
	}

	if(conv==0)
	{
		mbs->UO() <<"\n"<< "All the eigenvalues converged!"<<"\n";
		return 1;
	}
	mbs->UO() <<"\n"<< "NOT all eigenvalues converged!"<<"\n";
	return 0;

}