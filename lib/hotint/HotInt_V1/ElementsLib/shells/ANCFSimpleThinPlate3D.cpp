#include "ANCFSimpleThinPlate3D.h"

void ANCFSimpleThinPlate3D::SetANCFSimpleThinPlate3D(int body, const TArray<int>& nodes, int material, 
				double size1, double size2, double thickness, const Vector3D& color)
{
	this->SetFiniteElementGeneric(body, nodes, material, color);

	q_ref.SetLen(2*SOS());
	q_ref.FillWithZeros();
	for (int i = 1; i <= NNodes(); i++)
	{
		ANCFNodeS1S2_3D& n = (ANCFNodeS1S2_3D&)mbs->GetNode(nodes(i));

		// TODO: change, johannes doesn't want low level routines
		memcpy(this->q_ref.GetVecPtr()+0+(i-1)*DOFPerNode(), n.Pos().GetVecPtr(), 3*sizeof(double));
		memcpy(this->q_ref.GetVecPtr()+3+(i-1)*DOFPerNode(), n.GetRefSlope1().GetVecPtr() , 3*sizeof(double));
		memcpy(this->q_ref.GetVecPtr()+6+(i-1)*DOFPerNode(), n.GetRefSlope2().GetVecPtr() , 3*sizeof(double));
	}

	this->size1 = size1;
	this->size2 = size2;
	this->jac(1) = 2. / size1; this->jac(2) = 2. / size2;
	this->jac(3) = 4. / (size1*size1); this->jac(4) = 4. / (size1*size2); this->jac(5) = 4. / (size2*size2);
	this->thickness = thickness;
}

void ANCFSimpleThinPlate3D::ElementDefaultConstructorInitialization()
{
	elementname = GetElementSpec();
	nodes.SetLen(4); nodes.SetAll(0);
	size1 = 0.; size2 = 0.; thickness = 0.;
	x_init.SetLen(2*SOS()); x_init.FillWithZeros();
	materialnum = 0;
}

void ANCFSimpleThinPlate3D::EvalF2(Vector& f, double t) 
{
	//add loads, constraint forces and other things in parent function:
	Body3D::EvalF2(f,t);
	TMStartTimer(22);

	ConstVector<FEmaxDOF> fadd;
	fadd.SetLen(SOS());
	fadd.SetAll(0);

	if (GetGeometricNonlinearityStatus() == GNS_Linear)
	{
		EvalF2GeomLin(fadd, t);
	}
	else
	{
		EvalF2GeomNonlin(fadd, t);
	}

	f -= fadd;
}; 

void ANCFSimpleThinPlate3D::EvalF2GeomLin(Vector& f, double t) 
{
	// first cache xg for faster access
	SetXG_cached(TCD_reference_configuration);
	ConstVector<36> temp(SOS());
	Vector3D temp3d;

	double E = GetMaterial().YoungsModulus();
	double nu = GetMaterial().PoissonRatio();
	Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
				E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
				0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

	ConstVector<36> q; q.SetLen(SOS());
	for (int i = 1; i <= SOS(); i++)
	{
		q(i) = (XG_dc(i, TCD_compute) - XG_dc(i, TCD_reference_configuration));
	}

	for (IntegrationPointsIterator ip(this->integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		const Vector2D& p_loc = ip.Point2D();		
		Matrix grad;
		GetGrad(ip, grad);

		Vector3D rx = GetDPosDAlpha(p_loc, 1, grad, TCD_cached);
		Vector3D ry = GetDPosDAlpha(p_loc, 2, grad, TCD_cached);
		Vector3D rxx = GetDDPosDAlphaDBeta(p_loc, 1, grad, TCD_cached);
		Vector3D rxy = GetDDPosDAlphaDBeta(p_loc, 2, grad, TCD_cached);
		Vector3D ryy = GetDDPosDAlphaDBeta(p_loc, 3, grad, TCD_cached);

		Vector3D n = rx.Cross(ry);
		double norm = n.Norm();
		double norm2 = norm*norm;
		double norm3 = norm*norm2;
		n /= norm;

		ConstMatrix<36*3> delta_eps, delta_kappa;
		delta_eps.SetSize(SOS(), Dim());
		delta_kappa.SetSize(SOS(), Dim());

		double fac_eps = ip.Weight() * thickness * 0.25*size1*size2;
		double fac_kappa = fac_eps * thickness*thickness / 12.;

		for (int i = 1; i <= NS(); i++)
		{
			double dsdx = grad(i,1);
			double dsdy = grad(i,2);

			double ddsdxx  = grad(i,3);
			double ddsdxdy = grad(i,4);
			double ddsdyy  = grad(i,5);

			int ii = (i-1)*Dim();
			delta_eps(ii + 1, 1) = rx(1) * dsdx;
			delta_eps(ii + 2, 1) = rx(2) * dsdx;
			delta_eps(ii + 3, 1) = rx(3) * dsdx;

			delta_eps(ii + 1, 2) = ry(1) * dsdy;
			delta_eps(ii + 2, 2) = ry(2) * dsdy;
			delta_eps(ii + 3, 2) = ry(3) * dsdy;

			delta_eps(ii + 1, 3) = (ry(1) * dsdx + rx(1) * dsdy);
			delta_eps(ii + 2, 3) = (ry(2) * dsdx + rx(2) * dsdy);
			delta_eps(ii + 3, 3) = (ry(3) * dsdx + rx(3) * dsdy);
		}

		for (int i = 1; i <= SOS(); i++)
		{
			// The vectors drxdqi, drydqi, drxxdqi, drxydqi, dryydqi have a single non-zero component. for this reason, 
			// vector products are replaced by scalar operations. as a consequence, no vectors are needed any longer
			double dsdxi = grad((i-1)/3+1,1);		// drxdqi
			double dsdyi = grad((i-1)/3+1,2);		// drydqi
			double ddsdxxi  = grad((i-1)/3+1,3);	// drxxdqi
			double ddsdxdyi = grad((i-1)/3+1,4);	// drxydqi
			double ddsdyyi  = grad((i-1)/3+1,5);	// dryydqi

			//Vector3D drxdqi, drydqi;
			int ii = (i - 1) % 3 + 1;

			// Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);
			Vector3D dndqi;
			dndqi((ii  ) % 3 + 1) = -dsdyi * ry((ii+1) % 3 + 1) + dsdxi * rx((ii+1) % 3 + 1);
			dndqi((ii+1) % 3 + 1) =	 dsdyi * ry((ii  ) % 3 + 1) - dsdxi * rx((ii  ) % 3 + 1);
			
			delta_kappa(i,1) = -3./norm3 * (n*dndqi) * (n*rxx) + 1./norm3 * (dndqi*rxx) + 1./norm2 * (n(ii) * ddsdxxi);
			delta_kappa(i,2) = -3./norm3 * (n*dndqi) * (n*ryy) + 1./norm3 * (dndqi*ryy) + 1./norm2 * (n(ii) * ddsdyyi);
			delta_kappa(i,3) = -6./norm3 * (n*dndqi) * (n*rxy) + 2./norm3 * (dndqi*rxy) + 2./norm2 * (n(ii) * ddsdxdyi);
		}

		Vector3D epsC, kappaC; 
		MultTp(delta_eps, q, temp3d);
		Mult(C, temp3d, epsC);

		MultTp(delta_kappa, q, temp3d);
		Mult(C, temp3d, kappaC);

		Mult(delta_eps, epsC, temp);
		temp *= fac_eps;
		f += temp;

		Mult(delta_kappa, kappaC, temp);
		temp *= fac_kappa;
		f += temp;
	}

}; 

void ANCFSimpleThinPlate3D::EvalF2GeomNonlin(Vector& f, double t) 
{
#if 1

// first cache xg for faster access
	SetXG_cached();
	ConstVector<36> temp(SOS());

	double E = GetMaterial().YoungsModulus();
	double nu = GetMaterial().PoissonRatio();
	Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
				E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
				0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

	for (IntegrationPointsIterator ip(this->integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		const Vector2D& p_loc = ip.Point2D();		
		Matrix grad;
		GetGrad(ip, grad);

		Vector3D rx = GetDPosDAlpha(p_loc, 1, grad, TCD_cached);
		Vector3D ry = GetDPosDAlpha(p_loc, 2, grad, TCD_cached);
		Vector3D rxx = GetDDPosDAlphaDBeta(p_loc, 1, grad, TCD_cached);
		Vector3D rxy = GetDDPosDAlphaDBeta(p_loc, 2, grad, TCD_cached);
		Vector3D ryy = GetDDPosDAlphaDBeta(p_loc, 3, grad, TCD_cached);

		Vector3D eps, kappa;
		eps(1) = 0.5*(rx*rx - 1.);
		eps(2) = 0.5*(ry*ry - 1.);
		eps(3) = rx*ry;

		Vector3D n = rx.Cross(ry);
		double norm = n.Norm();
		double norm2 = norm*norm;
		double norm3 = norm*norm2;
		n.Normalize();

		kappa(1) = 1./norm2 * (n*rxx);
		kappa(2) = 1./norm2 * (n*ryy);
		kappa(3) = 1./norm2 * (n*rxy) * 2.;

		Vector3D epsC(	C(1,1)*eps(1)   + C(1,2)*eps(2),	C(2,1)*eps(1)   + C(2,2)*eps(2),	C(3,3)*eps(3));
		Vector3D kappaC(C(1,1)*kappa(1) + C(1,2)*kappa(2),	C(2,1)*kappa(1) + C(2,2)*kappa(2),	C(3,3)*kappa(3));

		ConstMatrix<36*3> delta_eps(SOS(), Dim()), delta_kappa(SOS(), Dim());

		double fac_eps = ip.Weight() * thickness * 0.25*size1*size2;
		double fac_kappa = fac_eps * thickness*thickness / 12.;

		for (int i = 1; i <= NS(); i++)
		{
			double dsdx = grad(i,1);
			double dsdy = grad(i,2);

			double ddsdxx  = grad(i,3);
			double ddsdxdy = grad(i,4);
			double ddsdyy  = grad(i,5);

			delta_eps((i-1)*Dim() + 1, 1) = rx(1) * dsdx;
			delta_eps((i-1)*Dim() + 2, 1) = rx(2) * dsdx;
			delta_eps((i-1)*Dim() + 3, 1) = rx(3) * dsdx;

			delta_eps((i-1)*Dim() + 1, 2) = ry(1) * dsdy;
			delta_eps((i-1)*Dim() + 2, 2) = ry(2) * dsdy;
			delta_eps((i-1)*Dim() + 3, 2) = ry(3) * dsdy;

			delta_eps((i-1)*Dim() + 1, 3) = (ry(1) * dsdx + rx(1) * dsdy);
			delta_eps((i-1)*Dim() + 2, 3) = (ry(2) * dsdx + rx(2) * dsdy);
			delta_eps((i-1)*Dim() + 3, 3) = (ry(3) * dsdx + rx(3) * dsdy);
		}

		for (int i = 1; i <= SOS(); i++)
		{
			double dsdx = grad((i-1)/3+1,1);
			double dsdy = grad((i-1)/3+1,2);
			double ddsdxx  = grad((i-1)/3+1,3);
			double ddsdxdy = grad((i-1)/3+1,4);
			double ddsdyy  = grad((i-1)/3+1,5);

			Vector3D drxdqi, drydqi, drxxdqi, drxydqi, dryydqi;
			int j = (i - 1) % 3 + 1;
			drxdqi(j) = dsdx;
			drydqi(j) = dsdy;

			drxxdqi(j) = ddsdxx;
			drxydqi(j) = ddsdxdy;
			dryydqi(j) = ddsdyy;

			Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);	
			
			double n_dndqi = n*dndqi;
			delta_kappa(i,1) = -3./norm3 * n_dndqi * (n*rxx) + 1./norm3 * (dndqi*rxx) + 1./norm2 * (n * drxxdqi);
			delta_kappa(i,2) = -3./norm3 * n_dndqi * (n*ryy) + 1./norm3 * (dndqi*ryy) + 1./norm2 * (n * dryydqi);
			delta_kappa(i,3) = -6./norm3 * n_dndqi * (n*rxy) + 2./norm3 * (dndqi*rxy) + 2./norm2 * (n * drxydqi);
		}

		Mult(delta_eps, epsC, temp);
		temp *= fac_eps;
		f += temp;

		Mult(delta_kappa, kappaC, temp);
		temp *= fac_kappa;
		f += temp;
	}
#else
	for (IntegrationPointsIterator ip(this->integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		const Vector2D& p_loc = ip.Point2D();

		double E = GetMaterial().YoungsModulus();
		double nu = GetMaterial().PoissonRatio();
		Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
					E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
					0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

		Vector3D epsC = C*GetMidPlaneStrains(p_loc, TCD_compute);
		Vector3D kappaC = C*GetCurvatures(p_loc, TCD_compute);

		Matrix delta_eps, delta_kappa;
		GetDeltaEps(p_loc, delta_eps, TCD_compute);
		GetDeltaKappa(p_loc, delta_kappa, TCD_compute);

		Vector temp(SOS());
		Mult(delta_eps, epsC, temp);
		f += temp * (ip.Weight() * thickness * 0.25*size1*size2);

		temp.FillWithZeros();
		Mult(delta_kappa, kappaC, temp);
		f += temp * (Cub(thickness) / 12. * ip.Weight() * 0.25*size1*size2);
	}
#endif
}; 

void ANCFSimpleThinPlate3D::EvalM(Matrix& m, double t)
{
	if (massmatrix.Getcols() == SOS())
	{
		m = massmatrix;
		return;
	}
	else
	{
		Matrix H(SOS(), Dim());

		for (IntegrationPointsIterator ip(this->integrationRuleMass); !ip.IsEnd(); ++ip)
		{
			const Vector2D& p_loc = ip.Point2D();

			for (int i = 1; i <= NS(); i++)
			{
				for (int j = 1; j <= Dim(); j++)
				{
					H(Dim()*(i-1)+j,j) = GetS0(p_loc, i);
				}
			}

			m += GetMaterial().Density() * ip.Weight() * thickness * (H*H.GetTp()) * 0.25*size1*size2;
		}
		massmatrix = m;
	}
}

void ANCFSimpleThinPlate3D::StiffnessMatrix(Matrix& m)
{
	if (GetGeometricNonlinearityStatus() == GNS_Linear)
	{
		StiffnessMatrixLin(m);
	}
	else
	{
		StiffnessMatrixNonlin(m);
	}
}

void ANCFSimpleThinPlate3D::StiffnessMatrixLin(Matrix& m)
{
	Vector3D temp3d;

	if (stiffnessmatrix.Getcols() == SOS())
	{
		//m.LinkWith(SOS(),SOS(), stiffnessmatrix.GetMatPtr());	// ist das gut sooo?
		m = stiffnessmatrix;
	}
	else
	{
		m.SetAll(0.);
		SetXG_cached(TCD_reference_configuration);

		double E = GetMaterial().YoungsModulus();
		double nu = GetMaterial().PoissonRatio();
		Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
					E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
					0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

		for (IntegrationPointsIterator ip(this->integrationRuleStiffness); !ip.IsEnd(); ++ip)
		{
			const Vector2D& p_loc = ip.Point2D();		
			Matrix grad;
			GetGrad(ip, grad);

			Vector3D rx = GetDPosDAlpha(p_loc, 1, grad, TCD_cached);
			Vector3D ry = GetDPosDAlpha(p_loc, 2, grad, TCD_cached);
			Vector3D rxx = GetDDPosDAlphaDBeta(p_loc, 1, grad, TCD_cached);
			Vector3D rxy = GetDDPosDAlphaDBeta(p_loc, 2, grad, TCD_cached);
			Vector3D ryy = GetDDPosDAlphaDBeta(p_loc, 3, grad, TCD_cached);

			Vector3D n = rx.Cross(ry);
			double norm = n.Norm();
			double norm2 = norm*norm;
			double norm3 = norm*norm2;
			n /= norm;

			double fac_eps = ip.Weight() * thickness * 0.25*size1*size2;
			double fac_kappa = fac_eps * thickness*thickness / 12.;
			
			for (int i = 1; i <= FlexDOF(); i++)
			{
				// The vectors drxdqi, drydqi, drxxdqi, drxydqi, dryydqi have a single non-zero component. for this reason, 
				// vector products are replaced by scalar operations. as a consequence, no vectors are needed any longer
				double dsdxi = grad((i-1)/3+1,1);		// drxdqi
				double dsdyi = grad((i-1)/3+1,2);		// drydqi
				double ddsdxxi  = grad((i-1)/3+1,3);	// drxxdqi
				double ddsdxdyi = grad((i-1)/3+1,4);	// drxydqi
				double ddsdyyi  = grad((i-1)/3+1,5);	// dryydqi

				//Vector3D drxdqi, drydqi;
				//Vector3D drxxdqi, drxydqi, dryydqi;
				int ii = (i - 1) % 3 + 1;
				Vector3D depsdqi, dkappadqi;
				depsdqi(1) = dsdxi * rx(ii);
				depsdqi(2) = dsdyi * ry(ii);
				depsdqi(3) = dsdxi * ry(ii) + rx(ii) * dsdyi;

				// Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);
				Vector3D dndqi;
				dndqi((ii  ) % 3 + 1) = -dsdyi * ry((ii+1) % 3 + 1) + dsdxi * rx((ii+1) % 3 + 1);
				dndqi((ii+1) % 3 + 1) =	 dsdyi * ry((ii  ) % 3 + 1) - dsdxi * rx((ii  ) % 3 + 1);
			
				double n_dndqi = (n*dndqi);
				double n_rxx = (n*rxx);
				double n_ryy = (n*ryy);
				double n_rxy = (n*rxy);

				dkappadqi(1) = -3./norm3 * n_dndqi * n_rxx + 1./norm3 * (dndqi*rxx) + 1./norm2 * (n(ii) * ddsdxxi);
				dkappadqi(2) = -3./norm3 * n_dndqi * n_ryy + 1./norm3 * (dndqi*ryy) + 1./norm2 * (n(ii) * ddsdyyi);
				dkappadqi(3) = -6./norm3 * n_dndqi * n_rxy + 2./norm3 * (dndqi*rxy) + 2./norm2 * (n(ii) * ddsdxdyi);

				for (int j = 1; j <= i; j++)
				{
					// The vectors drxdqj, drydqj, drxxdqj, drxydqj, dryydqj have a single non-zero component. for this reason, 
					// vector products are replaced by scalar operations. as a consequence, no vectors are needed any longer
					double dsdxj = grad((j-1)/3+1,1);		// drxdqj
					double dsdyj = grad((j-1)/3+1,2);		// drydqj
					double ddsdxxj  = grad((j-1)/3+1,3);	// drxxdqj
					double ddsdxdyj = grad((j-1)/3+1,4);	// drxydqj
					double ddsdyyj  = grad((j-1)/3+1,5);	// dryydqj

					int jj = (j - 1) % 3 + 1;
					Vector3D depsdqj, dkappadqj;
					depsdqj(1) = dsdxj * rx(jj);
					depsdqj(2) = dsdyj * ry(jj);
					depsdqj(3) = dsdxj * ry(jj) + rx(jj) * dsdyj;

					// Vector3D dndqj = drxdqj.Cross(ry) + rx.Cross(drydqj);
					Vector3D dndqj;
					dndqj((jj  ) % 3 + 1) = -dsdyj * ry((jj+1) % 3 + 1) + dsdxj * rx((jj+1) % 3 + 1);
					dndqj((jj+1) % 3 + 1) =	 dsdyj * ry((jj  ) % 3 + 1) - dsdxj * rx((jj  ) % 3 + 1);

					double n_dndqj = (n*dndqj);

					dkappadqj(1) = -3./norm3 * n_dndqj * n_rxx + 1./norm3 * (dndqj*rxx) + 1./norm2 * (n(jj) * ddsdxxj);
					dkappadqj(2) = -3./norm3 * n_dndqj * n_ryy + 1./norm3 * (dndqj*ryy) + 1./norm2 * (n(jj) * ddsdyyj);
					dkappadqj(3) = -6./norm3 * n_dndqj * n_rxy + 2./norm3 * (dndqj*rxy) + 2./norm2 * (n(jj) * ddsdxdyj);

					Mult(C, depsdqi, temp3d);
					m(i,j) -= (depsdqj*temp3d) * fac_eps;
					Mult(C, dkappadqi, temp3d);
					m(i,j) -= (dkappadqj*temp3d) * fac_kappa;
				}
			}
		}
		for (int i = 1; i <= FlexDOF(); i++)
		{
			for (int j = 1; j < i; j++)
			{
				m(j,i) = m(i,j);
			}
		};
	}
}

void ANCFSimpleThinPlate3D::StiffnessMatrixNonlin(Matrix& m)
{
	m.SetAll(0.);
	SetXG_cached();

	double E = GetMaterial().YoungsModulus();
	double nu = GetMaterial().PoissonRatio();
	Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
				E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
				0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

	for (IntegrationPointsIterator ip(this->integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		const Vector2D& p_loc = ip.Point2D();		
		ConstMatrix<12*5> grad;
		GetGrad(ip, grad);

		Vector3D rx = GetDPosDAlpha(p_loc, 1, grad, TCD_cached);
		Vector3D ry = GetDPosDAlpha(p_loc, 2, grad, TCD_cached);
		Vector3D rxx = GetDDPosDAlphaDBeta(p_loc, 1, grad, TCD_cached);
		Vector3D rxy = GetDDPosDAlphaDBeta(p_loc, 2, grad, TCD_cached);
		Vector3D ryy = GetDDPosDAlphaDBeta(p_loc, 3, grad, TCD_cached);

		Vector3D eps, kappa;
		eps(1) = 0.5*(rx*rx - 1.);
		eps(2) = 0.5*(ry*ry - 1.);
		eps(3) = rx*ry;

		Vector3D n = rx.Cross(ry);
		double norm = n.Norm();
		double norm2 = norm*norm;
		double norm3 = norm*norm2;
		double norm4 = norm2*norm2;
		n.Normalize();

		kappa(1) = (n*rxx);
		kappa(2) = (n*ryy);
		kappa(3) = (n*rxy) * 2.;
		kappa /= norm2;

		Vector3D epsC(	C(1,1)*eps(1)   + C(1,2)*eps(2),	C(2,1)*eps(1)   + C(2,2)*eps(2),	C(3,3)*eps(3));
		Vector3D kappaC(C(1,1)*kappa(1) + C(1,2)*kappa(2),	C(2,1)*kappa(1) + C(2,2)*kappa(2),	C(3,3)*kappa(3));

		double fac_eps = ip.Weight() * thickness * 0.25*size1*size2;
		double fac_kappa = fac_eps*thickness*thickness / 12.;
		
		for (int i = 1; i <= FlexDOF(); i++)
		{
			double dsdxi = grad((i-1)/3+1,1);
			double dsdyi = grad((i-1)/3+1,2);
			double ddsdxxi  = grad((i-1)/3+1,3);
			double ddsdxdyi = grad((i-1)/3+1,4);
			double ddsdyyi  = grad((i-1)/3+1,5);

			Vector3D drxdqi, drydqi;
			Vector3D drxxdqi, drxydqi, dryydqi;
			int ii = (i - 1) % 3 + 1;
			drxdqi(ii) = dsdxi;
			drydqi(ii) = dsdyi;
			drxxdqi(ii) = ddsdxxi;
			drxydqi(ii) = ddsdxdyi;
			dryydqi(ii) = ddsdyyi;

			Vector3D depsdqi, dkappadqi;
			depsdqi(1) = drxdqi * rx;
			depsdqi(2) = drydqi * ry;
			depsdqi(3) = drxdqi * ry + rx * drydqi;

			Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);

			double n_dndqi = (n*dndqi);
			double n_rxx = (n*rxx);
			double n_ryy = (n*ryy);
			double n_rxy = (n*rxy);

			dkappadqi(1) = -3./norm3 * n_dndqi * n_rxx + 1./norm3 * (dndqi*rxx) + 1./norm2 * (n * drxxdqi);
			dkappadqi(2) = -3./norm3 * n_dndqi * n_ryy + 1./norm3 * (dndqi*ryy) + 1./norm2 * (n * dryydqi);
			dkappadqi(3) = -6./norm3 * n_dndqi * n_rxy + 2./norm3 * (dndqi*rxy) + 2./norm2 * (n * drxydqi);

			Vector3D depsdqiC(	C(1,1)*depsdqi(1)   + C(1,2)*depsdqi(2),	C(2,1)*depsdqi(1)   + C(2,2)*depsdqi(2),	C(3,3)*depsdqi(3));
			Vector3D dkappadqiC(C(1,1)*dkappadqi(1) + C(1,2)*dkappadqi(2),	C(2,1)*dkappadqi(1) + C(2,2)*dkappadqi(2),	C(3,3)*dkappadqi(3));

			for (int j = 1; j <= i; j++)
			{
				double dsdxj = grad((j-1)/3+1,1);
				double dsdyj = grad((j-1)/3+1,2);
				double ddsdxxj  = grad((j-1)/3+1,3);
				double ddsdxdyj = grad((j-1)/3+1,4);
				double ddsdyyj  = grad((j-1)/3+1,5);

				Vector3D drxdqj, drydqj;
				Vector3D drxxdqj, drxydqj, dryydqj;
				int jj = (j - 1) % 3 + 1;
				drxdqj(jj) = dsdxj;
				drydqj(jj) = dsdyj;
				drxxdqj(jj) = ddsdxxj;
				drxydqj(jj) = ddsdxdyj;
				dryydqj(jj) = ddsdyyj;
				
				Vector3D ddeps, ddkappa;
				ddeps(1) = drxdqi*drxdqj;
				ddeps(2) = drydqi*drydqj;
				ddeps(3) = drxdqi*drydqj + drxdqj*drydqi;

				Vector3D depsdqj, dkappadqj;
				depsdqj(1) = drxdqj * rx;
				depsdqj(2) = drydqj * ry;
				depsdqj(3) = drxdqj * ry + rx * drydqj;

				Vector3D dndqj = drxdqj.Cross(ry) + rx.Cross(drydqj);
				Vector3D ddndqidqj = drxdqi.Cross(drydqj) + drxdqj.Cross(drydqi);
				double dndqi_dndqj = (dndqj*dndqi);
				double n_ddndqidqj = (n*ddndqidqj);
				double n_dndqj = (n*dndqj);

				dkappadqj(1) = -3./norm3 * n_dndqj * n_rxx + 1./norm3 * (dndqj*rxx) + 1./norm2 * (n * drxxdqj);
				dkappadqj(2) = -3./norm3 * n_dndqj * n_ryy + 1./norm3 * (dndqj*ryy) + 1./norm2 * (n * dryydqj);
				dkappadqj(3) = -6./norm3 * n_dndqj * n_rxy + 2./norm3 * (dndqj*rxy) + 2./norm2 * (n * drxydqj);

				
				ddkappa(1) = 1./norm3 * ((ddndqidqj*rxx) + (dndqi*drxxdqj)
					- 3./norm * dndqi_dndqj*n_rxx + 15./norm * n_dndqj*n_dndqi*n_rxx - 3. * n_ddndqidqj*n_rxx
					- 3./norm * n_dndqi*(dndqj*rxx) - 3. * n_dndqi*(n*drxxdqj) 
					+ 1. * (dndqj*drxxdqi) - 3. * n_dndqj*(n*drxxdqi));
				
				ddkappa(2) = 1./norm3 * ((ddndqidqj*ryy) + 1. * (dndqi*dryydqj)
					- 3./norm * dndqi_dndqj*n_ryy + 15./norm * n_dndqj*n_dndqi*n_ryy - 3. * n_ddndqidqj*n_ryy
					- 3./norm * n_dndqi*(dndqj*ryy) - 3. * n_dndqi*(n*dryydqj) 
					+ 1. * (dndqj*dryydqi) - 3. * n_dndqj*(n*dryydqi));
				
				ddkappa(3) = (1./norm3 * ((ddndqidqj*rxy) + 1. * (dndqi*drxydqj)
					- 3./norm * dndqi_dndqj*n_rxy + 15./norm * n_dndqj*n_dndqi*n_rxy - 3. * n_ddndqidqj*n_rxy
					- 3./norm * n_dndqi*(dndqj*rxy) - 3. * n_dndqi*(n*drxydqj) 
					+ 1. * (dndqj*drxydqi) - 3. * n_dndqj*(n*drxydqi))) * 2.;

				m(i,j) -= (ddeps*epsC + depsdqj*depsdqiC) * fac_eps;
				m(i,j) -= (ddkappa*kappaC + dkappadqj*dkappadqiC) * fac_kappa;
			}
		}
	}
	for (int i = 1; i <= FlexDOF(); i++)
	{
		for (int j = 1; j < i; j++)
		{
			m(j,i) = m(i,j);
		}
	};
}


double ANCFSimpleThinPlate3D::GetJacInvDS(const Vector3D& p_loc, Matrix& jacinvDS) const
{
	// a bit misused: store shape function derivatives here
	jacinvDS.SetSize(NS(), 5);

	for (int i = 1; i <= NS(); i++)
	{
		jacinvDS(i,1) = GetDS0((Vector2D&)p_loc,i,1) * 2. / size1;
		jacinvDS(i,2) = GetDS0((Vector2D&)p_loc,i,2) * 2. / size2;
		jacinvDS(i,3) = GetDDS0((Vector2D&)p_loc,i,1) * 4. / (size1*size1);
		jacinvDS(i,4) = GetDDS0((Vector2D&)p_loc,i,2) * 4. / (size1*size2);
		jacinvDS(i,5) = GetDDS0((Vector2D&)p_loc,i,3) * 4. / (size2*size2);
	}
	return size1*size2 / 4.;
}

Matrix3D ANCFSimpleThinPlate3D::GetElementJacobian(const Vector2D& p_loc) const
{
	/*Vector3D refposxi = GetDPosDAlpha(p_loc, 1, TCD_reference_configuration);
	Vector3D refposeta = GetDPosDAlpha(p_loc, 2, TCD_reference_configuration);
	Vector3D n = GetNormal(p_loc, TCD_reference_configuration);
	Matrix3D jac(refposxi.X(), refposeta.X(), n.X(),
			refposxi.Y(), refposeta.Y(), n.Y(),
			refposxi.Z(), refposeta.Z(), n.Z());

	refposxi.Normalize();
	refposxi.GramSchmidt(refposeta);
	refposeta.Normalize();

	Matrix3D rot(refposxi.X(), refposxi.Y(), refposxi.Z(),
			refposeta.X(), refposeta.Y(), refposeta.Z(),
			n.X(), n.Y(), n.Z());

	return rot*jac;*/
	return Matrix3D(size1/2.,	0.,			0.,
					0.,			size2/2.,	0., 
					0.,			0.,			1.);
}

Vector3D ANCFSimpleThinPlate3D::GetMidPlaneStrains(const Vector2D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D e;
	Vector3D rx = GetPosx(p_loc, flag);
	Vector3D ry = GetPosy(p_loc, flag);

	e(1) = 0.5*(rx*rx - 1.);
	e(2) = 0.5*(ry*ry - 1.);
	e(3) = rx*ry;
	return e;
}

Vector3D ANCFSimpleThinPlate3D::GetCurvatures(const Vector2D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D k;
	Vector3D rxx = GetPosxx(p_loc, flag);
	Vector3D rxy = GetPosxy(p_loc, flag);
	Vector3D ryy = GetPosyy(p_loc, flag);

	Vector3D rx = GetPosx(p_loc, flag);
	Vector3D ry = GetPosy(p_loc, flag);
	
	Vector3D n = rx.Cross(ry);
	double d = pow(n*n, 1.5);

	k(1) = 1./d * (n*rxx);
	k(2) = 1./d * (n*ryy);
	k(3) = 1./d * (n*rxy) * 2.;
	return k;
}

Vector3D ANCFSimpleThinPlate3D::GetStrains(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D v = GetMidPlaneStrains((Vector2D&)p_loc, flag);
	if (p_loc.Z() != 0)
	{
		v += GetCurvatures((Vector2D&)p_loc, flag) * (p_loc.Z() * thickness / 2.); 
	}
	return v;
}

Vector3D ANCFSimpleThinPlate3D::GetStresses(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	double E = GetMaterial().YoungsModulus();
	double nu = GetMaterial().PoissonRatio();
	Matrix3D C(	E/(1.-nu*nu),		E*nu/(1.-nu*nu),	0.,
				E*nu/(1.-nu*nu),	E/(1.-nu*nu),		0.,
				0.,					0.,					E*(1.-nu)/(2.*(1.-nu*nu)) );

	Vector3D strains = GetStrains(p_loc, flag);
	
	Matrix2D eps, stress; 
	StrainVectorToMatrix2D(eps, strains);
	//GetMaterial().ComputeStressFromStrain2D(eps, stress, Vector(0));
	return C*strains;
}

void ANCFSimpleThinPlate3D::GetDeltaEps(const Vector2D& p_loc, Matrix& delta_eps, TComputeDrawInitFlag flag) const
{
	delta_eps.SetSize(SOS(), Dim());
	delta_eps.FillWithZeros();

	Vector3D rx = GetPosx(p_loc, flag);
	Vector3D ry = GetPosy(p_loc, flag);

	for (int i = 1; i <= NS(); i++)
	{
		double s0x =  GetDS0(p_loc, i, 1) * 2. / size1;
		double s0y =  GetDS0(p_loc, i, 2) * 2. / size2;

		delta_eps((i-1)*Dim() + 1, 1) = rx(1) * s0x;
		delta_eps((i-1)*Dim() + 2, 1) = rx(2) * s0x;
		delta_eps((i-1)*Dim() + 3, 1) = rx(3) * s0x;

		delta_eps((i-1)*Dim() + 1, 2) = ry(1) * s0y;
		delta_eps((i-1)*Dim() + 2, 2) = ry(2) * s0y;
		delta_eps((i-1)*Dim() + 3, 2) = ry(3) * s0y;

		delta_eps((i-1)*Dim() + 1, 3) = (ry(1) * s0x + rx(1) * s0y);
		delta_eps((i-1)*Dim() + 2, 3) = (ry(2) * s0x + rx(2) * s0y);
		delta_eps((i-1)*Dim() + 3, 3) = (ry(3) * s0x + rx(3) * s0y);
	}
}

void ANCFSimpleThinPlate3D::GetDeltaKappa(const Vector2D& p_loc, Matrix& delta_kappa, TComputeDrawInitFlag flag) const
{
	delta_kappa.SetSize(SOS(), Dim());
	delta_kappa.FillWithZeros();

	Vector3D rx = GetPosx(p_loc, flag);
	Vector3D ry = GetPosy(p_loc, flag);

	Vector3D rxx = GetPosxx(p_loc, flag);
	Vector3D rxy = GetPosxy(p_loc, flag);
	Vector3D ryy = GetPosyy(p_loc, flag);

	Vector3D n = rx.Cross(ry);
	double n2 = n*n;

	Vector3D ex(1.,0.,0.), ey(0.,1.,0.), ez(0.,0.,1.);

	for (int i = 1; i <= NS(); i++)
	{
		double s0x = GetDS0(p_loc, i, 1) * 2. / size1;
		double s0y = GetDS0(p_loc, i, 2) * 2. / size2;

		double s0xx = GetDDS0(p_loc, i, 1) * 4. / (size1*size1);
		double s0xy = GetDDS0(p_loc, i, 2) * 4. / (size1*size2);
		double s0yy = GetDDS0(p_loc, i, 3) * 4. / (size2*size2);

		delta_kappa((i-1)*Dim() + 1, 1) += 1./pow(n2, 3./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * rxx;
		delta_kappa((i-1)*Dim() + 2, 1) += 1./pow(n2, 3./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * rxx;
		delta_kappa((i-1)*Dim() + 3, 1) += 1./pow(n2, 3./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * rxx;

		delta_kappa((i-1)*Dim() + 1, 1) -= 3./pow(n2, 5./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * n * (n * rxx);
		delta_kappa((i-1)*Dim() + 2, 1) -= 3./pow(n2, 5./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * n * (n * rxx);
		delta_kappa((i-1)*Dim() + 3, 1) -= 3./pow(n2, 5./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * n * (n * rxx);

		delta_kappa((i-1)*Dim() + 1, 1) += 1./pow(n2, 3./2.) * n.X() * s0xx;
		delta_kappa((i-1)*Dim() + 2, 1) += 1./pow(n2, 3./2.) * n.Y() * s0xx;
		delta_kappa((i-1)*Dim() + 3, 1) += 1./pow(n2, 3./2.) * n.Z() * s0xx;


		delta_kappa((i-1)*Dim() + 1, 2) += 1./pow(n2, 3./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * ryy;
		delta_kappa((i-1)*Dim() + 2, 2) += 1./pow(n2, 3./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * ryy;
		delta_kappa((i-1)*Dim() + 3, 2) += 1./pow(n2, 3./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * ryy;

		delta_kappa((i-1)*Dim() + 1, 2) -= 3./pow(n2, 5./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * n * (n * ryy);
		delta_kappa((i-1)*Dim() + 2, 2) -= 3./pow(n2, 5./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * n * (n * ryy);
		delta_kappa((i-1)*Dim() + 3, 2) -= 3./pow(n2, 5./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * n * (n * ryy);

		delta_kappa((i-1)*Dim() + 1, 2) += 1./pow(n2, 3./2.) * n.X() * s0yy;
		delta_kappa((i-1)*Dim() + 2, 2) += 1./pow(n2, 3./2.) * n.Y() * s0yy;
		delta_kappa((i-1)*Dim() + 3, 2) += 1./pow(n2, 3./2.) * n.Z() * s0yy;


		delta_kappa((i-1)*Dim() + 1, 3) += 2./pow(n2, 3./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * rxy;
		delta_kappa((i-1)*Dim() + 2, 3) += 2./pow(n2, 3./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * rxy;
		delta_kappa((i-1)*Dim() + 3, 3) += 2./pow(n2, 3./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * rxy;

		delta_kappa((i-1)*Dim() + 1, 3) -= 6./pow(n2, 5./2.) * (s0x * ex.Cross(ry) + s0y * rx.Cross(ex)) * n * (n * rxy);
		delta_kappa((i-1)*Dim() + 2, 3) -= 6./pow(n2, 5./2.) * (s0x * ey.Cross(ry) + s0y * rx.Cross(ey)) * n * (n * rxy);
		delta_kappa((i-1)*Dim() + 3, 3) -= 6./pow(n2, 5./2.) * (s0x * ez.Cross(ry) + s0y * rx.Cross(ez)) * n * (n * rxy);

		delta_kappa((i-1)*Dim() + 1, 3) += 2./pow(n2, 3./2.) * n.X() * s0xy;
		delta_kappa((i-1)*Dim() + 2, 3) += 2./pow(n2, 3./2.) * n.Y() * s0xy;
		delta_kappa((i-1)*Dim() + 3, 3) += 2./pow(n2, 3./2.) * n.Z() * s0xy;
	}
}

Vector3D ANCFSimpleThinPlate3D::GetPos_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D p;

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0((const Vector2D &)p_loc, j);
		p(1) += s0 * XG_dc((j-1)*Dim()+1, flag);
		p(2) += s0 * XG_dc((j-1)*Dim()+2, flag);
		p(3) += s0 * XG_dc((j-1)*Dim()+3, flag);
	}
	if (p_loc.Z() != 0)
	{
		// double thickness = GetThicknessAtPoint((const Vector2D &)ploc, xg.IsDrawConfiguration());
		p += GetNormal((const Vector2D&)p_loc, flag) * (p_loc.Z() * thickness / 2.);
	}
	
	return p;
}

Vector3D ANCFSimpleThinPlate3D::GetDisplacement_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	return GetPos_dc(p_loc, flag) - GetRefPos(p_loc);
}

Vector3D ANCFSimpleThinPlate3D::GetVel_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D v0;

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0((const Vector2D&)p_loc, j);
		v0(1) += s0 * XG_dc((j-1)*Dim()+1+SOS(), flag);
		v0(2) += s0 * XG_dc((j-1)*Dim()+2+SOS(), flag);
		v0(3) += s0 * XG_dc((j-1)*Dim()+3+SOS(), flag);
	}

	Vector3D dndt;
	for (int i = 1; i <= SOS(); i++)
	{
		dndt += GetDNormalDqi((const Vector2D&)p_loc, i, flag) * XG_dc(i+SOS(), flag);
	}
	
	return v0 + dndt * (p_loc.Z() * thickness / 2.);
}

Vector3D ANCFSimpleThinPlate3D::GetNormal(const Vector2D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D n = GetPosx(p_loc, flag).Cross(GetPosy(p_loc, flag));
	n.Normalize();
	return n;
}

Vector3D ANCFSimpleThinPlate3D::GetDNormalDqi(const Vector2D& p_loc, int dof, TComputeDrawInitFlag flag) const
{
	Vector3D dndqi;

	Vector3D rx = GetPosx(p_loc, flag);
	Vector3D ry = GetPosy(p_loc, flag);
	
	Vector3D n = rx.Cross(ry);
	double norm = n.Norm();
	n.Normalize();

	double s0x = GetDS0(p_loc, (dof-1)/3+1, 1) * 2. / size1;
	double s0y = GetDS0(p_loc, (dof-1)/3+1, 2) * 2. / size2;

	Vector3D drxdqi, drydqi;
	drxdqi((dof - 1) % 3 + 1) = s0x;
	drydqi((dof - 1) % 3 + 1) = s0y;

	dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);
	dndqi -= n * (dndqi * n);

	return dndqi * (1./norm);
}

Vector3D ANCFSimpleThinPlate3D::GetDPosDAlpha(const Vector2D& p_loc, int alpha, TComputeDrawInitFlag flag) const
{
	Vector3D dp;
	for (int j = 1; j <= NS(); j++)
	{
		double ds0 = GetDS0(p_loc,j,alpha) * jac(alpha);
		dp(1) += ds0 * XG_dc((j-1)*Dim()+1, flag);
		dp(2) += ds0 * XG_dc((j-1)*Dim()+2, flag);
		dp(3) += ds0 * XG_dc((j-1)*Dim()+3, flag);
	}
	return dp;
}

Vector3D ANCFSimpleThinPlate3D::GetDDPosDAlphaDBeta(const Vector2D& p_loc, int alpha_beta, TComputeDrawInitFlag flag) const
{
	Vector3D ddp;
	for (int j = 1; j <= NS(); j++)
	{
		double dds0 = GetDDS0(p_loc,j,alpha_beta) * jac(2+alpha_beta);
		ddp(1) += dds0 * XG_dc((j-1)*Dim()+1, flag);
		ddp(2) += dds0 * XG_dc((j-1)*Dim()+2, flag);
		ddp(3) += dds0 * XG_dc((j-1)*Dim()+3, flag);
	}
	return ddp;
}

Vector3D ANCFSimpleThinPlate3D::GetDPosPDAlpha(const Vector2D& p_loc, int alpha, TComputeDrawInitFlag flag) const
{
	Vector3D dpP;
	for (int j = 1; j <= NS(); j++)
	{
		double ds0 = GetDS0(p_loc,j,alpha) * jac(alpha);
		dpP(1) += ds0 * XGP_dc((j-1)*Dim()+1, flag);
		dpP(2) += ds0 * XGP_dc((j-1)*Dim()+2, flag);
		dpP(3) += ds0 * XGP_dc((j-1)*Dim()+3, flag);
	}
	return dpP;
}

Vector3D ANCFSimpleThinPlate3D::GetDPosDAlpha(const Vector2D& p_loc, int alpha, Matrix& grad, TComputeDrawInitFlag flag) const
{
	Vector3D dp;
	for (int j = 1; j <= NS(); j++)
	{
		double ds0 = grad(j,alpha);
		dp(1) += ds0 * XG_dc((j-1)*Dim()+1, flag);
		dp(2) += ds0 * XG_dc((j-1)*Dim()+2, flag);
		dp(3) += ds0 * XG_dc((j-1)*Dim()+3, flag);
	}
	return dp;
}

Vector3D ANCFSimpleThinPlate3D::GetDDPosDAlphaDBeta(const Vector2D& p_loc, int alpha_beta, Matrix& grad, TComputeDrawInitFlag flag) const
{
	Vector3D ddp;
	for (int j = 1; j <= NS(); j++)
	{
		double dds0 = grad(j,2+alpha_beta);
		ddp(1) += dds0 * XG_dc((j-1)*Dim()+1, flag);
		ddp(2) += dds0 * XG_dc((j-1)*Dim()+2, flag);
		ddp(3) += dds0 * XG_dc((j-1)*Dim()+3, flag);
	}
	return ddp;
}

void ANCFSimpleThinPlate3D::GetdPosdqT(const Vector3D& p_loc, Matrix& m)
{
	m.SetSize(SOS(), Dim());
	m.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0((Vector2D&)p_loc, i);
		m((i-1)*Dim() + 1, 1) = s;
		m((i-1)*Dim() + 2, 2) = s;
		m((i-1)*Dim() + 3, 3) = s;
	}
	// TODO: out of plane position!!!!
}

Matrix3D ANCFSimpleThinPlate3D::GetRotMatrix_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D rx = GetPosx((Vector2D&)p_loc, flag);
	Vector3D ry = GetPosy((Vector2D&)p_loc, flag);

	// get an orthogonal frame in the deformed configuration, that almost(!) coincides with the actual rotated frame 
	// in case the deformation is not too large
	Vector3D v3 = rx.Cross(ry);
	Vector3D v2 = v3.Cross(rx);

	// normalize the vectors
	rx.Normalize(); v2.Normalize(); v3.Normalize();

	Matrix3D rot;
	rot.Set(rx, v2, v3);

	return rot;
}

Matrix3D ANCFSimpleThinPlate3D::GetRotMatrixP_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const
{
	Vector3D rx = GetPosx((Vector2D&)p_loc, flag);
	Vector3D ry = GetPosy((Vector2D&)p_loc, flag);

	Vector3D rxP = GetPosxP((Vector2D&)p_loc, flag);
	Vector3D ryP = GetPosyP((Vector2D&)p_loc, flag);

	Vector3D& v1 = rx;
	Vector3D v3 = rx.Cross(ry);
	Vector3D v2 = v3.Cross(rx);

	double v1_norm = v1.Norm();
	double v2_norm = v2.Norm();
	double v3_norm = v3.Norm();

	Vector3D v1P = (1./v1_norm) * rxP;
	v1P -= rx * ((rxP*rx) / (v1_norm*v1_norm)); 

	Vector3D nP = rxP.Cross(ry) + rx.Cross(ryP);
	Vector3D v3P = (1./v3_norm) * nP;
	v3P -= v3 * ((nP*v3) / (v3_norm*v3_norm));

	Vector3D v2P = v1P.Cross(v3) + v1.Cross(v3P);

	Matrix3D rotP;
	rotP.Set(v1P, v2P, v3P);

	return rotP;
}

void ANCFSimpleThinPlate3D::GetdRotdqT(const Vector3D& p_loc, Matrix& m)
{
	m.SetSize(SOS(), Dim());
	m.FillWithZeros();

	Vector3D rx = GetPosx((Vector2D&)p_loc, TCD_compute);
	Vector3D ry = GetPosy((Vector2D&)p_loc, TCD_compute);

	// get an orthogonal frame in the deformed configuration, that almost(!) coincides with the actual rotated frame 
	// in case the deformation is not too large
	Vector3D& v1 = rx;
	Vector3D v3 = rx.Cross(ry);
	Vector3D v2 = v3.Cross(rx);

	double v1_norm = v1.Norm();
	double v2_norm = v2.Norm();
	double v3_norm = v3.Norm();

	for (int i = 1; i <= SOS(); i++)
	{
		double dsdxi = GetDS0((Vector2D&)p_loc, (i-1)/3+1, 1) * jac(1);
		double dsdyi = GetDS0((Vector2D&)p_loc, (i-1)/3+1, 2) * jac(2);

		Vector3D drxdqi, drydqi;
		int ii = (i - 1) % 3 + 1;
		drxdqi(ii) = dsdxi;
		drydqi(ii) = dsdyi;

		Vector3D dv1dqi = (1./v1_norm) * drxdqi;
		dv1dqi -= rx * ((dv1dqi*rx) / (v1_norm*v1_norm)); 

		Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);
		Vector3D dv3dqi = (1./v3_norm) * dndqi;
		dv3dqi -= v3 * ((dv3dqi*v3) / (v3_norm*v3_norm));

		Vector3D dv2dqi = (1./v2_norm) * (dndqi.Cross(v1) + v3.Cross(drxdqi));
		dv2dqi -= v2 * ((dv2dqi*v2) / (v2_norm*v2_norm));

		Vector3D dthetadi 
			= (v1*(dv2dqi*v3))*(1./(v1_norm*v3_norm)) 
			+ (v2*(dv3dqi*v1))*(1./(v1_norm*v2_norm)) 
			+ (v3*(dv1dqi*v2))*(1./(v2_norm*v3_norm));

		m(i,1) = dthetadi(1);
		m(i,2) = dthetadi(2);
		m(i,3) = dthetadi(3);
	}
}

void ANCFSimpleThinPlate3D::GetdRotvdqT(const Vector3D& v_loc, const Vector3D& p_loc, Matrix& m)
{
	m.SetSize(SOS(), Dim());
	m.FillWithZeros();

	Vector3D rx = GetPosx((Vector2D&)p_loc, TCD_compute);
	Vector3D ry = GetPosy((Vector2D&)p_loc, TCD_compute);

	// get an orthogonal frame in the deformed configuration, that almost(!) coincides with the actual rotated frame 
	// in case the deformation is not too large
	Vector3D& v1 = rx;
	Vector3D v3 = rx.Cross(ry);
	Vector3D v2 = v3.Cross(rx);

	double v1_norm = v1.Norm();
	double v2_norm = v2.Norm();
	double v3_norm = v3.Norm();

	for (int i = 1; i <= SOS(); i++)
	{
		double dsdxi = GetDS0((Vector2D&)p_loc, (i-1)/3+1, 1) * jac(1);
		double dsdyi = GetDS0((Vector2D&)p_loc, (i-1)/3+1, 2) * jac(2);

		Vector3D drxdqi, drydqi;
		int ii = (i - 1) % 3 + 1;
		drxdqi(ii) = dsdxi;
		drydqi(ii) = dsdyi;

		Vector3D dv1dqi = (1./v1_norm) * drxdqi;
		dv1dqi -= rx * ((dv1dqi*rx) / (v1_norm*v1_norm)); 

		Vector3D dndqi = drxdqi.Cross(ry) + rx.Cross(drydqi);
		Vector3D dv3dqi = (1./v3_norm) * dndqi;
		dv3dqi -= v3 * ((dv3dqi*v3) / (v3_norm*v3_norm));

		Vector3D dv2dqi = (1./v2_norm) * (dndqi.Cross(v1) + v3.Cross(drxdqi));
		dv2dqi -= v2 * ((dv2dqi*v2) / (v2_norm*v2_norm));

		//v1 /= v1_norm; v2 /= v2_norm; v3 /= v3_norm;

		m(i,1) = dv1dqi(1) * v_loc(1) + dv2dqi(1) * v_loc(2) + dv3dqi(1) * v_loc(3);
		m(i,2) = dv1dqi(2) * v_loc(1) + dv2dqi(2) * v_loc(2) + dv3dqi(2) * v_loc(3);
		m(i,3) = dv1dqi(3) * v_loc(1) + dv2dqi(3) * v_loc(2) + dv3dqi(3) * v_loc(3);
	}
}


Vector3D ANCFSimpleThinPlate3D::GetNodeLocPos(int local_node_number) const
{
	switch(local_node_number)
	{
	case 1: return Vector3D(-1.,-1.,0.);
	case 2: return Vector3D( 1.,-1.,0.);
	case 3: return Vector3D( 1., 1.,0.);
	case 4: return Vector3D(-1., 1.,0.);
	default: assert(0); return Vector3D();
	}
}

void ANCFSimpleThinPlate3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor>& variables)
{
	FiniteElementGeneric<Body3D>::GetAvailableFieldVariables(variables);

	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,	FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress, FieldVariableDescriptor::FVCI_y, true);
}

double ANCFSimpleThinPlate3D::GetFieldVariableValue_dc(const FieldVariableDescriptor& fvd, const Vector3D& p_loc, TComputeDrawInitFlag flag)
{
	switch (fvd.VariableType())
	{
		case FieldVariableDescriptor::FVT_position: return fvd.GetComponent(GetPos_dc(p_loc, flag));
		case FieldVariableDescriptor::FVT_displacement: return fvd.GetComponent(GetDisplacement_dc(p_loc, flag));
		case FieldVariableDescriptor::FVT_velocity: return fvd.GetComponent(GetVel_dc(p_loc, flag));
		case FieldVariableDescriptor::FVT_total_strain: 
			{
				Vector3D eps = GetStrains(p_loc, flag);
				if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_x && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_x)
				{
					return eps(1);
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_y)
				{
					return eps(2);
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_x && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_y)
				{
					return eps(3)/2.;
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_magnitude)
				{
					return 0.;
				}
			}
		case FieldVariableDescriptor::FVT_stress: 
			{
				Vector3D sigma = GetStresses(p_loc, flag);
				if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_x && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_x)
				{
					return sigma(1);
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_y)
				{
					return sigma(2);
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_x && fvd.ComponentIndex2() == FieldVariableDescriptor::FVCI_y)
				{
					return sigma(3);
				}
				else if (fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_magnitude)
				{
					return sqrt(sigma(1)*sigma(1)+sigma(2)*sigma(2)-sigma(1)*sigma(2)+3.*sigma(3)*sigma(3));
				}
			}
	}
	return FIELD_VARIABLE_NO_VALUE;
}

// shape functions as they are, accessible per index nsf = 1 .. 12
// local coordinates q1, q2 vary from -1 to +1
double ANCFSimpleThinPlate3D::GetS0(const Vector2D& p_loc, int nsf) const
{
	double xi = p_loc.X();
	double eta = p_loc.Y();

	switch(nsf)
	{
	case 1: return xi*eta/2.0-xi*xi*xi*eta/8.0-xi*eta*eta*eta/8.0+1.0/4.0+xi*xi*xi/8.0+eta*eta
		*eta/8.0-3.0/8.0*eta-3.0/8.0*xi;
	case 2: return size1*xi*xi*eta/16.0+size1/16.0+size1*xi*eta/16.0-size1*xi/16.0-size1*xi*xi/16.0-size1*
		xi*xi*xi*eta/16.0+size1*xi*xi*xi/16.0-size1*eta/16.0;
	case 3: return size2/16.0-size2*xi*eta*eta*eta/16.0-size2*eta/16.0-size2*eta*eta/16.0+size2*eta*eta*eta/
		16.0-size2*xi/16.0+size2*xi*eta/16.0+size2*xi*eta*eta/16.0;
	case 4: return -xi*eta/2.0+xi*xi*xi*eta/8.0+xi*eta*eta*eta/8.0+1.0/4.0-xi*xi*xi/8.0+eta*
		eta*eta/8.0-3.0/8.0*eta+3.0/8.0*xi;
	case 5: return -size1/16.0-size1*xi/16.0+size1*xi*xi/16.0+size1*xi*xi*xi/16.0+size1*eta/16.0+size1*xi*
		eta/16.0-size1*xi*xi*eta/16.0-size1*xi*xi*xi*eta/16.0;
	case 6: return size2/16.0-size2*eta/16.0+size2*xi/16.0-size2*xi*eta/16.0-size2*eta*eta/16.0-size2*xi*eta*
		eta/16.0+size2*eta*eta*eta/16.0+size2*xi*eta*eta*eta/16.0;
	case 7: return xi*eta/2.0-xi*xi*xi*eta/8.0-xi*eta*eta*eta/8.0+1.0/4.0-xi*xi*xi/8.0-eta*eta
		*eta/8.0+3.0/8.0*eta+3.0/8.0*xi;
	case 8: return -size1/16.0-size1*eta/16.0-size1*xi/16.0-size1*xi*eta/16.0+size1*xi*xi/16.0+size1*xi*xi*
		eta/16.0+size1*xi*xi*xi/16.0+size1*xi*xi*xi*eta/16.0;
	case 9: return -size2/16.0-size2*eta/16.0+size2*eta*eta/16.0-size2*xi/16.0-size2*xi*eta/16.0+size2*xi*eta*
		eta/16.0+size2*eta*eta*eta/16.0+size2*xi*eta*eta*eta/16.0;
	case 10: return -xi*eta/2.0+xi*xi*xi*eta/8.0+xi*eta*eta*eta/8.0+1.0/4.0+xi*xi*xi/8.0-eta*
		eta*eta/8.0+3.0/8.0*eta-3.0/8.0*xi;
	case 11: return size1/16.0+size1*eta/16.0-size1*xi/16.0-size1*xi*eta/16.0-size1*xi*xi/16.0-size1*xi*xi*
		eta/16.0+size1*xi*xi*xi/16.0+size1*xi*xi*xi*eta/16.0;
	case 12: return -size2/16.0-size2*eta/16.0+size2*eta*eta/16.0+size2*eta*eta*eta/16.0+size2*xi/16.0+size2*xi*
			eta/16.0-size2*xi*eta*eta/16.0-size2*xi*eta*eta*eta/16.0;
	}
	assert(0);
	return 0;
}

// derivatives of the shape functions with respect to the local coordinates; alpha = 1,2
double ANCFSimpleThinPlate3D::GetDS0(const Vector2D& p_loc, int nsf, int alpha) const
{
	double xi = p_loc.X();
	double eta = p_loc.Y();

	if(alpha == 1)
	{
		// derivatives with respect to xi
		switch(nsf)
		{
		case 1: return eta/2.0-3.0/8.0*xi*xi*eta-eta*eta*eta/8.0+3.0/8.0*xi*xi-3.0/8.0;
		case 2: return size1*xi*eta/8.0+size1*eta/16.0-size1/16.0-size1*xi/8.0-3.0/16.0*size1*xi*xi*eta+3.0/16.0*size1*xi*xi;
		case 3: return -size2*eta*eta*eta/16.0-size2/16.0+size2*eta/16.0+size2*eta*eta/16.0;
		case 4: return -eta/2.0+3.0/8.0*xi*xi*eta+eta*eta*eta/8.0-3.0/8.0*xi*xi+3.0/8.0;
		case 5: return -size1/16.0+size1*xi/8.0+3.0/16.0*size1*xi*xi+size1*eta/16.0-size1*xi*eta/8.0-3.0/16.0*size1*xi*xi*eta;
		case 6: return size2/16.0-size2*eta/16.0-size2*eta*eta/16.0+size2*eta*eta*eta/16.0;
		case 7: return eta/2.0-3.0/8.0*xi*xi*eta-eta*eta*eta/8.0-3.0/8.0*xi*xi+3.0/8.0;
		case 8: return -size1/16.0-size1*eta/16.0+size1*xi/8.0+size1*xi*eta/8.0+3.0/16.0*size1*xi*xi+3.0/16.0*size1*xi*xi*eta;
		case 9: return -size2/16.0-size2*eta/16.0+size2*eta*eta/16.0+size2*eta*eta*eta/16.0;
		case 10: return -eta/2.0+3.0/8.0*xi*xi*eta+eta*eta*eta/8.0+3.0/8.0*xi*xi-3.0/8.0;
		case 11: return -size1/16.0-size1*eta/16.0-size1*xi/8.0-size1*xi*eta/8.0+3.0/16.0*size1*xi*xi+3.0/16.0*size1*xi*xi*eta;
		case 12: return size2/16.0+size2*eta/16.0-size2*eta*eta/16.0-size2*eta*eta*eta/16.0;
		}
	}
	else if(alpha == 2)
	{
		// derivatives with respect to eta
		switch(nsf)
		{
		case 1: return xi/2.0-xi*xi*xi/8.0-3.0/8.0*xi*eta*eta+3.0/8.0*eta*eta-3.0/8.0;
		case 2: return size1*xi*xi/16.0+size1*xi/16.0-size1*xi*xi*xi/16.0-size1/16.0;
		case 3: return -3.0/16.0*size2*xi*eta*eta-size2/16.0-size2*eta/8.0+3.0/16.0*size2*eta*eta+size2*xi/16.0+size2*xi*eta/8.0;
		case 4: return -xi/2.0+xi*xi*xi/8.0+3.0/8.0*xi*eta*eta+3.0/8.0*eta*eta-3.0/8.0;
		case 5: return size1/16.0+size1*xi/16.0-size1*xi*xi/16.0-size1*xi*xi*xi/16.0;
		case 6: return -size2/16.0-size2*xi/16.0-size2*eta/8.0-size2*xi*eta/8.0+3.0/16.0*size2*eta*eta+3.0/16.0*size2*xi*eta*eta;
		case 7: return xi/2.0-xi*xi*xi/8.0-3.0/8.0*xi*eta*eta-3.0/8.0*eta*eta+3.0/8.0;
		case 8: return -size1/16.0-size1*xi/16.0+size1*xi*xi/16.0+size1*xi*xi*xi/16.0;
		case 9: return -size2/16.0+size2*eta/8.0-size2*xi/16.0+size2*xi*eta/8.0+3.0/16.0*size2*eta*eta+3.0/16.0*size2*xi*eta*eta;
		case 10: return -xi/2.0+xi*xi*xi/8.0+3.0/8.0*xi*eta*eta-3.0/8.0*eta*eta+3.0/8.0;
		case 11: return size1/16.0-size1*xi/16.0-size1*xi*xi/16.0+size1*xi*xi*xi/16.0;
		case 12: return -size2/16.0+size2*eta/8.0+3.0/16.0*size2*eta*eta+size2*xi/16.0-size2*xi*eta/8.0-3.0/16.0*size2*xi*eta*eta;
		}
	}
	assert(0);
	return 0;
}

// second order derivatives of the shape functions with respect to the local coordinates
// alpha_beta = 1 (11), 2 (12), 3 (22)
double ANCFSimpleThinPlate3D::GetDDS0(const Vector2D& p_loc, int nsf, int alpha_beta) const
{
	double xi = p_loc.X();
	double eta = p_loc.Y();

	if(alpha_beta == 1)
	{
	// derivatives d^2 S / d xi^2
		switch(nsf)
		{
		case 1: return -3.0/4.0*xi*eta+3.0/4.0*xi;
		case 2: return -3.0/8.0*size1*xi*eta+size1*eta/8.0+3.0/8.0*size1*xi-size1/8.0;
		case 3: return 0.0;
		case 4: return 3.0/4.0*xi*eta-3.0/4.0*xi;
		case 5: return size1/8.0+3.0/8.0*size1*xi-size1*eta/8.0-3.0/8.0*size1*xi*eta;
		case 6: return 0.0;
		case 7: return -3.0/4.0*xi*eta-3.0/4.0*xi;
		case 8: return size1/8.0+size1*eta/8.0+3.0/8.0*size1*xi+3.0/8.0*size1*xi*eta;
		case 9: return 0.0;
		case 10: return 3.0/4.0*xi*eta+3.0/4.0*xi;
		case 11: return -size1/8.0-size1*eta/8.0+3.0/8.0*size1*xi+3.0/8.0*size1*xi*eta;
		case 12: return 0.0;
		}
	}
	else if(alpha_beta == 2)
	{
	// derivatives d^2 S / (d xi d eta)
		switch(nsf)
		{
		case 1: return 1.0/2.0-3.0/8.0*xi*xi-3.0/8.0*eta*eta;
		case 2: return -3.0/16.0*size1*xi*xi+size1*xi/8.0+size1/16.0;
		case 3: return size2/16.0-3.0/16.0*size2*eta*eta+size2*eta/8.0;
		case 4: return -1.0/2.0+3.0/8.0*xi*xi+3.0/8.0*eta*eta;
		case 5: return size1/16.0-size1*xi/8.0-3.0/16.0*size1*xi*xi;
		case 6: return -size2/16.0-size2*eta/8.0+3.0/16.0*size2*eta*eta;
		case 7: return 1.0/2.0-3.0/8.0*xi*xi-3.0/8.0*eta*eta;
		case 8: return -size1/16.0+size1*xi/8.0+3.0/16.0*size1*xi*xi;
		case 9: return -size2/16.0+size2*eta/8.0+3.0/16.0*size2*eta*eta;
		case 10: return -1.0/2.0+3.0/8.0*xi*xi+3.0/8.0*eta*eta;
		case 11: return -size1/16.0-size1*xi/8.0+3.0/16.0*size1*xi*xi;
		case 12: return size2/16.0-size2*eta/8.0-3.0/16.0*size2*eta*eta;
		}
	}
	else if(alpha_beta == 3)
	{
	// derivatives d^2 S / d eta^2
		switch(nsf)
		{
		case 1: return -3.0/4.0*xi*eta+3.0/4.0*eta;
		case 2: return 0.0;
		case 3: return -3.0/8.0*size2*xi*eta+size2*xi/8.0-size2/8.0+3.0/8.0*size2*eta;
		case 4: return 3.0/4.0*xi*eta+3.0/4.0*eta;
		case 5: return 0.0;
		case 6: return -size2/8.0-size2*xi/8.0+3.0/8.0*size2*eta+3.0/8.0*size2*xi*eta;
		case 7: return -3.0/4.0*xi*eta-3.0/4.0*eta;
		case 8: return 0.0;
		case 9: return size2/8.0+size2*xi/8.0+3.0/8.0*size2*eta+3.0/8.0*size2*xi*eta;
		case 10: return 3.0/4.0*xi*eta-3.0/4.0*eta;
		case 11: return 0.0;
		case 12: return size2/8.0+3.0/8.0*size2*eta-size2*xi/8.0-3.0/8.0*size2*xi*eta;
		}
	}
	assert(0);
	return 0;
}

void ANCFSimpleThinPlate3D::DefineIntegrationRule(IntegrationRule& integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_ThinPlate);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
	{
		if(integrationRule.settings.geometricNonlinearityStatus == GNS_Linear)
			ruleOrder = 4;	// previously 6
		else
			ruleOrder = 8;
	}
	else if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
	{
		if(integrationRule.settings.geometricNonlinearityStatus == GNS_Linear)
			ruleOrder = 2;
		else
			ruleOrder = 4;
	}
	else
		ruleOrder = 4;
	assert(ruleOrder != 0);

	IntegrationRule::DefineIntegrationRuleSquare(integrationRule, ruleOrder);
}

void ANCFSimpleThinPlate3D::AddSurfacePressure(Vector& f, double pressure, int dir)
{
	// in the linear case, the outer normal always points in the same (i.e., initial) direction
	if (this->GetGeometricNonlinearityStatus() == GNS_Linear)
	{
		this->SetXG_cached(TCD_reference_configuration);
	}
	else
	{
		this->SetXG_cached(TCD_compute);
	}

	Matrix H(SOS(), Dim());

	if (dir != 4 && dir != 5) 
	{
		assert(0);
	}
	if (dir == 4) 
	{	
		pressure *= -1.;
	}

	for (IntegrationPointsIterator ip(this->integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		const Vector2D& p_loc = ip.Point2D();

		Vector3D n = GetNormal(p_loc, TCD_cached);

		for (int i = 1; i <= NS(); i++)
		{
			double fact = GetS0(p_loc, i) * pressure * ip.Weight() * 0.25*size1*size2;

			for (int j = 1; j <= Dim(); j++)
			{
				f(Dim()*(i-1)+j) += fact*n(j);
			}
		}
	}
}

Vector3D ANCFSimpleThinPlate3D::GetSurfaceNormalD(int dir)
{
	// always return normal to mid-surface since the plate is supposed to be thin
	int sign = (dir == 4 ? -1 : 1);
	if (dir == 4 || dir == 5)
	{
		Vector3D n = GetNormal(Vector2D(), TCD_draw_magnified);
		n *= sign;
		return n;
	}
	return Vector3D();
}

Vector3D ANCFSimpleThinPlate3D::GetDOFPosD(int idof) const
{
	int nodeStartDOF = ((idof-1)/9) * 9 + 1;

	return Vector3D(
		XG_dc(nodeStartDOF + 0, TCD_draw_magnified),
		XG_dc(nodeStartDOF + 1, TCD_draw_magnified),
		XG_dc(nodeStartDOF + 2, TCD_draw_magnified)
		);
}

Vector3D ANCFSimpleThinPlate3D::GetDOFDirD(int idof) const
{
	int dir = (idof-1) % 3;

	if (dir == 0) return Vector3D(1.,0.,0.);
	else if (dir == 1) return Vector3D(0.,1.,0.);
	else return Vector3D(0.,0.,1.);
}

void ANCFSimpleThinPlate3D::DrawElement()
{
	mbs->SetColor(col);
	double lx1 = GetMBS()->GetDOption(106);		// Shrinking factor
	double ly1 = GetMBS()->GetDOption(106);		// Shrinking factor
	double lz1 = GetMBS()->GetDOption(106)*GetMBS()->GetMagnifyYZ();
	int drawflat = GetMBS()->GetIOption(117);	//draw only element midplane

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	} 
	else if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}
	else if (!GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 3;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double tilex = GetMBS()->GetIOption(137);	// axis resolution
	if (!colormode)
		tilex = GetMBS()->GetIOption(136);		// axis tiling for faces

	double tiley = tilex;

	static TArray<Vector3D> points;
	points.SetLen((int)(tilex+1)*(int)(tiley+1));
	static TArray<Vector3D> normals;
	normals.SetLen((int)(tilex+1)*(int)(tiley+1));
	static TArray<double> vals;
	vals.SetLen((int)(tilex+1)*(int)(tiley+1));

	int starts = 1;
	int ends = 6;
	if (drawflat)
	{
		starts = 3;
		ends = 3;
	}

	for (int side = starts; side <= ends; side++)
	{
		points.SetLen(0);
		vals.SetLen(0);
		normals.SetLen(0);
		Vector3D p0, vx, vy;
		int tileyn = (int)tiley;
		if (side <= 2 || side >= 5)
			tileyn = GetMBS()->GetIOption(138);		// tiling in transverse direction

		switch(side)
		{
		case 1:
			{ //bottom
				p0 = Vector3D(-lx1,-ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		case 2:
			{ //top
				p0 = Vector3D(-lx1, ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,-2.*lz1/tileyn);
				break;
			}
		case 3:
			{ //front
				if (drawflat)
					p0 = Vector3D(-lx1,-ly1, 0);
				else
					p0 = Vector3D(-lx1,-ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,2.*ly1/tileyn,0);
				
				break;
			}
		case 4:
			{ //back
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,-2.*ly1/tileyn,0);
				break;
			}
		case 5:
			{ //left
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(0,-2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		case 6:
			{ //right
				p0 = Vector3D( lx1,-ly1,-lz1);
				vx = Vector3D(0,2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		}

		for (double iy = 0; iy <= tileyn+1e-10; iy++)
		{
			for (double ix = 0; ix <= tilex+1e-10; ix++)
			{
				Vector3D ploc = (p0+ix*vx+iy*vy);
				Vector3D v = GetPos_dc(ploc, TCD_draw_magnified);
				points.Add(GetPos_dc(ploc, TCD_draw_magnified));

				if (side == 3)
					normals.Add(GetNormal(ploc, TCD_draw_magnified));
				if (side == 4)
					normals.Add(-1.*GetNormal(ploc, TCD_draw_magnified));

				if (colormode)
				{
					vals.Add(GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, TCD_draw));
				}
				else
					vals.Add(0);
			}
		}
		for (int i=1; i <= normals.Length(); i++)
			points.Add(normals(i));

		// FitPointsInDrawingRegion(points, (int)tilex+1, (int)tileyn+1);		// in case the derived class wishes to control the drawing region

		mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);
	}
}