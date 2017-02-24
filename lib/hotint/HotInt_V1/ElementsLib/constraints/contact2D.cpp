//#**************************************************************
//#
//# filename:             contact2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17. October 2006
//# description:          contact formulations
//#                       
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
 

#include "element.h"
#include "constraint.h"
#include "rigid2d.h"
#include "geomelements.h"
#include "contact2d.h"
#include "femathhelperfunctions.h"

// AP 2013-01-11 had to add this function, since it is not found in MBSKernelLib
int Find(int i, const IVector& v)
{
	for (int j=1; j <= v.Length(); j++)
	{
		if (v(j) == i) return j;
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GeneralContact2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GeneralContact2D::BuildSearchTrees()
{
	TMStartTimer(24); //build search trees
	//if (testcnt >= 1) //original: >5 
	//{
	//	testcnt = 0;
		//+++++++++++++++++++++++++++++
		//master segments:
		Box2D box;
		box.Clear();
		for (int i=1; i <= NMC(); i++)
		{
			Box2D pbox;
			GeomElement* ge = GeomElem(i);
			pbox = ge->GetBoundingBox2D();
			//UO() << "tp1=" << ge->GetTPoint2D(1) << "\n";
			//UO() << "STpbox=" << pbox.PMin() << "-" << pbox.PMax() << "\n";
			box.Add(pbox);
		}
		//UO() << "STbox=" << box.PMin() << "-" << box.PMax() << "\n";
		box.Increase(bordersize);

		mastersearchtree.ResetSearchTree(searchtreeix, searchtreeiy, box); //define number of cubes and maximum size of search tree
		for (int i=1; i <= NMC(); i++)
		{
			Box2D pbox = GeomElem(i)->GetBoundingBox2D();

			pbox.Increase(bordersize);
			mastersearchtree.AddItem(pbox,i);
		}
	//}
	//else testcnt ++;

	TMStopTimer(24); //build search trees
}

const Body2D& GeneralContact2D::GetSlaveBody(int i) const 
{
	return (const Body2D&)mbs->GetElement(slaveelements(i));
}

Body2D& GeneralContact2D::GetSlaveBody(int i) 
{
	return (Body2D&)mbs->GetElement(slaveelements(i));
}

const Body2D& GeneralContact2D::GetMasterBody(int i) const
{
	//$ YV 2012-12-13: changed GeomElem(i)->GetBody2D() to GeomElem(i)->GetElement()
	return (Body2D&)(GeomElem(i)->GetElement());
}

Body2D& GeneralContact2D::GetMasterBody(int i) 
{
	return (Body2D&)(GeomElem(i)->GetElement());
}

Vector2D GeneralContact2D::GetSlaveNodePos(int i) const 
{
	if (slaveNODEmode == 0)
	{
		return GetSlaveBody(i).GetPos2D(loccoords(i));
	}
	else
	{
		return GetSlaveBody(i).GetNodePos2D(locnodenums(i));
	}
}

Vector2D GeneralContact2D::GetSlaveNodePosD(int i) const 
{
	if (slaveNODEmode == 0)
	{
		return GetSlaveBody(i).GetPos2DD(loccoords(i));
	}
	else
	{
		return GetSlaveBody(i).GetNodePos2DD(locnodenums(i));
	}
}

Vector2D GeneralContact2D::GetSlaveNodeVel(int i) const 
{
	if (slaveNODEmode == 0)
	{
		if(!GetMBS()->DoStaticComputation())
		{
			return GetSlaveBody(i).GetVel2D(loccoords(i));
		}
		else
		{
			double dt = GetMBS()->GetStepSize();

			double* oldptr = GetMBS()->GetXact().GetVecPtr();
			int oldlen = GetMBS()->GetXact().Length();
			
			const Vector &lsv = GetMBS()->GetLastSolVector();
			//GetMBS()->SetActState(lsv);
			mbs->SetActState(lsv);
			//GetMBS()->UO() << "Oldconfig = " << lsv << "\n";

			
			Vector2D p0 = GetSlaveBody(i).GetPos2D(loccoords(i));
			//GetMBS()->UO() << "p0 = " << p0 << "\n";


			//GetMBS()->SetActState(oldptr, oldlen);
			mbs->SetActState(oldptr, oldlen);

			Vector2D p1 = GetSlaveBody(i).GetPos2D(loccoords(i));
			//GetMBS()->UO() << "p1 = " << p1 << "\n";

			if (dt == 0) 
			{
				dt = 1;
				GetMBS()->UO() << "warning: getslavenodevel: dt=0\n";
			}
			Vector2D pret((p1(1)-p0(1))/dt,(p1(2)-p0(2))/dt);
			return pret;
		}
	}
	else
	{
		//***
		return GetSlaveBody(i).GetNodeVel2D(locnodenums(i));
	}
}

int GeneralContact2D::GetMasterIsActiveCondition(double t, int masterelem) const
{
	int cond = GetMasterIsActiveMathFunction(masterisactive(masterelem))->CheckCondition(GetMBS()->GetTime());
	if (t > 0.1 && cond != 0) 
	{
		GetMBS()->UO() << "warning, cond=" << cond << "\n";
	}

	return cond;
	//mathfunctions(masterisactive(j))->CheckCondition(GetMBS()->GetTime())
}

Vector2D GeneralContact2D::GetMasterLocPos(const Vector2D& pp, int j, int locind) const
{
	return GeomElem(j)->GetLocPos(locind, pp);
}

//approximate velocity at segment point pglob due linear interpolation
Vector2D GeneralContact2D::GetMasterSegVel(const Vector2D& pglob, int j, int locind) const 
{
	if(!GetMBS()->DoStaticComputation())
	{
		return GeomElem(j)->GetVel2D(GeomElem(j)->GetLocPos(locind, pglob));
	}
	else
	{
		double dt = GetMBS()->GetStepSize();

		double* oldptr = GetMBS()->GetXact().GetVecPtr();
		int oldlen = GetMBS()->GetXact().Length();
		//GetMBS()->SetActState(GetMBS()->GetLastSolVector());
		mbs->SetActState(GetMBS()->GetLastSolVector());

		Vector2D p0 = GeomElem(j)->GetLocPos(locind, pglob);

		//GetMBS()->SetActState(oldptr, oldlen);
		mbs->SetActState(oldptr, oldlen);

		Vector2D p1 = GeomElem(j)->GetLocPos(locind, pglob);
		if (dt == 0) 
		{
			dt = 1;
			GetMBS()->UO() << "warning: GetMasterSegVel: dt=0\n";
		}
		//return (p1-p0)/dt;
		Vector2D pret((p1(1)-p0(1))/dt,(p1(2)-p0(2))/dt);
		return pret;
	}
	
}

//get mindist and projected point pp
void GeneralContact2D::GetNearestMasterSegment(int slavenum, const Vector2D& p, int& masternum, int& locind, double& dist, Vector2D& pp)
{
	Box2D box;
	box.Add(p);
	//box.Increase(pointradius(slavenum));

	mastersearchtree.GetItemsInBox(box, itemlist);
	Vector2D p1, p2, pf;

	double mindist = 1e20;
	double d;
	int fseg = 0;
	locind = 0; 
	int flocind = 0;
	double gap;
	
	static IVector itemlist2;
	itemlist2.SetLen(0);
	//delete double elements in itemlist:
	for (int i=1; i <= itemlist.Length(); i++)
	{
		if (tempitemlist(itemlist(i)) == 0) 
		{
			itemlist2.Add(itemlist(i));
			tempitemlist(itemlist(i)) = 1;
		}
	}

/*
	for (int i = 1; i <= itemlist.Length(); i++)
	{
		int j = itemlist(i);
		*/
	for (int i = 1; i <= itemlist2.Length(); i++)
	{
		int j = itemlist2(i);
		tempitemlist(j) = 0;

		if (masterbodyindex(j) < slavebodyindex(slavenum) &&
				(masterisactive(j) == 0 || GetMasterIsActiveCondition(GetMBS()->GetTime(), j)) )
		{
			gap = GeomElem(j)->GetNearestPoint(p, locind, pf);//-pointradius(slavenum);
			d = Dist(p,pf);//-pointradius(slavenum);

			if (fabs(d) < fabs(mindist) 
				&& fabs(gap) < contactmaxdist) 
			{
				mindist = d;
				fseg = j;
				pp = pf;
				flocind = locind;
			}
		}
	}
	//dist = mindist;
	if (fseg)
	{
		dist = GeomElem(fseg)->GetNearestPoint(p, locind, pf);//-pointradius(slavenum);
	} 
	else dist = mindist;

	masternum = fseg;
}

void GeneralContact2D::GetNearestMasterSegments(int slavenum, const Vector2D& p, double circlerad, TArray<int>& a_masternum, 
		TArray<int>& a_locind, TArray<double>& a_dist, TArray<Vector2D>& a_pp)
{
	Box2D box;
	box.Add(p);
	box.Increase(circlerad);

	mastersearchtree.GetItemsInBox(box, itemlist);
	Vector2D p1, p2, pf;

	double d;
	int flocind = 0;
	double gap;

	int	locind = 0; 
	int masternum = 0;
	double dist = 1e20;
	Vector2D pp;

	//check for all elements which penetrate -> store
	//if no penetration, return nearest element!

	a_masternum.SetLen(0);
	a_locind.SetLen(0);
	a_dist.SetLen(0);
	a_pp.SetLen(0);

	
	static IVector itemlist2;
	itemlist2.SetLen(0);
	//delete double elements in itemlist:
	for (int i=1; i <= itemlist.Length(); i++)
	{
		if (tempitemlist(itemlist(i)) == 0) 
		{
			itemlist2.Add(itemlist(i));
			tempitemlist(itemlist(i)) = 1;
		}
	}

/*
	for (int i = 1; i <= itemlist.Length(); i++)
	{
		int j = itemlist(i);
		*/
	for (int i = 1; i <= itemlist2.Length(); i++)
	{
		int j = itemlist2(i);
		tempitemlist(j) = 0;
//PG BEGIN: oneway for SPH
		if (masterbodyindex(j) < slavebodyindex(slavenum) && (masterisactive(j) == 0 || GetMasterIsActiveCondition(GetMBS()->GetTime(), j)) )		
//PG END
		//if ((1 /*IsSPHContact(masterpodyindex(j), slavebodyindex(slavenum))*/ || masterbodyindex(j) < slavebodyindex(slavenum)) &&
		//		(masterisactive(j) == 0 || GetMasterIsActiveCondition(GetMBS()->GetTime(), j)) )
		{
			gap = GeomElem(j)->GetNearestPoint(p, locind, pf) - circlerad; //gap value??
			d = Dist(p,pf) - circlerad;

			if (fabs(d) < fabs(dist) && fabs(d) < contactmaxdist) 
			{
				dist = d;
				masternum = j;
				pp = pf;
				flocind = locind;
			}

			if (circlerad != 0 && d < 0	&& fabs(d) < contactmaxdist)  
			{
				//if (!Find(j, a_masternum)) 
				{
					int f = 0;
					int k = 1;
					while (k <= a_masternum.Length() && !f) 
					{
						//if (masterbodyindex(a_masternum(k)) == masterbodyindex(j)) f = k; 
						if (Dist2(pf, a_pp(k)) < 1e-16)
						{
							f = k; 
						}
						k++;
					}

					if (!f)
					{
						a_masternum.Add(j);
						a_dist.Add(d);
						a_pp.Add(pf);
						a_locind.Add(locind);
					}
					else
					{
						if (fabs(d) < fabs(a_dist(f)))
						{
							a_masternum(f) = j;
							a_dist(f) = d;
							a_pp(f) = pf;
							a_locind(f) = locind;
						}

					}
				}
			}
		}

	}
	if (masternum && a_masternum.Length() == 0 && circlerad == 0)
	{
		dist = GeomElem(masternum)->GetNearestPoint(p, locind, pf) - circlerad;
		//dist = Dist(p,pf) - circlerad;

		if (dist < 0)
		{
			a_dist.Add(dist);
			a_masternum.Add(masternum);
			a_pp.Add(pp);
			a_locind.Add(flocind);
		}
	}
	//otherwise leave arrays empty
}

void GeneralContact2D::GetTangentStickForce(int i, const Vector2D& t, double tvel, const Vector2D& pp, const Vector2D& lastpp, Vector2D& tforce)
{
  double friction_penalty = 0.1*c_contact(i);

	tforce = 0;
	if (islaststick(i) && isstick(i)) tforce = friction_penalty*(lastpp-pp);
	//else tforce = (-friction_penalty*tvel)*t;
}

double GeneralContact2D::GetContactForce(int i, double gap, double gapp)
{
	//double e = 0.95; //coeff of restitution; 0.95 ususally, with ambrosio
	double e = restitution_coeff; //coeff of restitution; 0.95 ususally
	double nh = hertzian_contact_param; //coeff for Hertzian contact; 1 works good

	double gi = gapp_init(i);
	if (fabs(gapp) > gapp_init(i)) gi = fabs(gapp);
	if (fabs(gi) < 1e-16) gi = contactmaxdist;

	switch (contactmode)
	{
	case 0: //Hertian contact
		return  -Sgn(gap)*c_contact(i)*(1.+0.75*(1.-Sqr(e))*gapp/gi)*pow(fabs(gap),nh); //XG(i) is compression force, Hertzian contact model
		break;
	case 1: //tied, linear elastic
		return  -c_contact(i)*gap;
		break;
	case 3: //linear contact model
		return  -c_contact(i)*gap;
		break;
	default: UO() << "ERROR: contactmode in GetContactForce\n"; return 0;
	}
}

void GeneralContact2D::PostprocessingStep() 
{
	//TMStartTimer(24);
	for (int i = 1; i <= NSC(); i++)
	{
		if (pointradius(i) == 0)
		{
			if (!iscontact(i)) gapp_init(i) = gapp_init_min(i);

			if (isfriction)
			{
				if (iscontact(i) && isstick(i)) 
				{
					//UO() << "stick" << i << "\n";
					if (islaststick(i) != isstick(i))
					{
						int j = isstick(i);

						Vector2D pglob = GetSlaveNodePos(i); 

						//store local mastersegment position:
						laststickp(i) = GetMasterLocPos(pglob, j, masterlocind(j));
						islaststick(i) = isstick(i); //master segment number
					}
				} 
				else
				{
					if (!iscontact(i))
					{
						islaststick(i) = 0;
					} 
					else
					{
						int j = iscontact(i);

						Vector2D pglob = GetSlaveNodePos(i);

						//store local mastersegment position:
						laststickp(i) = GetMasterLocPos(pglob, j, masterlocind(j));
					}

				}
				if (abs(slipdir(i)) == 2) slipdir(i) /= 2;
			}
		}
	}

	BuildSearchTrees();
	nlstepcnt = 0;
	//TMStopTimer(24);
}


double GeneralContact2D::PostNewtonStep(double t) 
{
	TMStartTimer(25);

	if (contactmode == 1 || contactmode == 2) //always contact
	{
		UO() << "ERROR: contact mode 1 and 2 not yet implemented for General Contact\n";
		/*
		rv = 0;
		for (i = 1; i <= NSC(); i++)
		{
		if (!iscontact(i))
		{
		iscontact(i) = 1;
		rv = 1;
		}
		}
		return rv;
		*/
	}

	// AP 18.12.2012: Use Memberfunction MaxNLIt() instead of int maxnlit, returns 3
	//int maxnlit = 3;
	if (nlstepcnt > /*maxnlit*/MaxNLIt()) 
	{
		TMStopTimer(25);
		return 0;
	}
	double nlerror = 0;

	for (int i = 1; i <= NMC(); i++)
	{
		mastertoslave[i-1].SetLen(0);
	}

	for (int i = 1; i <= NSC(); i++)
	{
		Vector2D p, pp;
		double gap;
		int mnum, mlocind;

		p = GetSlaveNodePos(i);

		int n=1;
		if (pointradius(i) == 0)
		{
			// get nearest segment or circle
			GetNearestMasterSegment(i, p, mnum, mlocind, gap, pp); //gap < 0 means penetration
			if (mnum == 0) n = 0;
		}
		else
		{
			GetNearestMasterSegments(i, p, pointradius(i), tempmasternum, templocind, tempdist, temp_pp);
			n=tempmasternum.Length();
		}

		//if (tempmasternum.Length() && t > 2.83)
		//	UO() << "master=" << tempmasternum << ", gaps=" << tempdist << "\n";


		for (int fi=1; fi <= n; fi++)
		{
		//if (mnum != 0)
		//{
			if (pointradius(i) != 0)
			{
				mnum = tempmasternum(fi);
				mlocind = templocind(fi);
				gap = tempdist(fi);
				pp = temp_pp(fi);
			}

			//UO() << "iscontact(" << i << ")=" << mnum << "\n";

			masterlocind(mnum) = mlocind;

			Vector2D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
			Vector2D mpmid = GeomElem(mnum)->GetPos2D(ploc);
			Vector2D n;
			if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
			else {n = p - mpmid; n.Normalize();}
			Vector2D tt(-n.Y(), n.X());

			double gapp;
			Vector2D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
			gapp = -((GetSlaveNodeVel(i) - mv) * n);      //project velocity to contact surface normal
			double tvel = ((GetSlaveNodeVel(i) - mv) * tt); //project velocity tangential to surface

			// AP this is not correct e.g. for GeomLine-master
			// AP instead do nothing, gap is computed correctly in GetNearestMasterSegment above!!!
			// gap = (p-mpmid)*n - pointradius(i);


			double c_force;
			if (islagrange) 
				c_force = XG(i);
			else
				c_force = GetContactForce(i, gap, gapp);

			if (pointradius(i) != 0)
			{
				if (Find(mnum,*multimastercontact(i)))
				{
					iscontact(i) = mnum;
				}
				else
				{
					iscontact(i) = 0;
					nlerror += 1;
				}
			}
			else
				if (iscontact(i) && mnum != iscontact(i)) 
				{ 
					iscontact(i) = mnum;
					nlerror += 1;
				}

			if (gap < 0 && iscontact(i) == 0 && nlstepcnt < MaxNLIt() )
			{
				iscontact(i) = mnum; //master element number
				nlerror += fabs(gap);
			}
			else 	if (gap > master_adhesion_coefficient(mnum) && iscontact(i) != 0)
			{
				iscontact(i) = 0;
				nlerror += fabs(c_force)/c_contact(i);
			}

			if (gapp_init(i) < fabs(gapp) && iscontact(i) != 0)
			{
				nlerror += fabs(gapp_init(i) - fabs(gapp))/Maximum(gapp_init(i),fabs(gapp));
				gapp_init(i) = fabs(gapp);
			}

			if (isfriction && pointradius(i) == 0)
			{
				if (isstick(i))
				{
					Vector2D t_force;
					Vector2D lastpp = GeomElem(isstick(i))->GetPos2D(laststickp(i));
					GetTangentStickForce(i, tt, tvel, pp, lastpp, t_force);

					//if (t >= 0.16 && t <= 0.1601) UO() << "slave" << i << ": gap=" << gap << ", cforce=" << c_force << ", tforce=" << t_force << "\n";

					
					if (t_force.Norm() > frictioncoeff_stick*c_force)
					{
						isstick(i) = 0;
						nlerror += fabs(t_force.Norm() - frictioncoeff_stick*c_force)/c_contact(i);
						if (abs(slipdir(i)) == 2) slipdir(i) = -slipdir(i) / 2;
					}
				}
				else
				{
					//isstick(i) = mnum;
					
					if ((double)slipdir(i)*tvel < 0)
					{
						//assumed slip direction not equal tangential slip velocity -> try stick or other velocity:
						isstick(i) = mnum;
						nlerror += fabs(tvel);
						slipdir(i) *= 2;
					}
				}
			}

			if (iscontact(i) != 0)
			{
				mastertoslave[mnum-1].Add(i);
			}

			if (pointradius(i) != 0)
			{
				tempmasternum(fi) = iscontact(i);
			}

			/*
			if (gap < 0.0 && (t >= 0.0))
			{
				UO() << "t" << t << ", gap" << i << "=" << gap << ", cforce=" << c_force << ", mnum=" << mnum;
				UO() << ", iscontact=" << iscontact(i) << ", err=" << nlerror << "\n";
				//UO() << "cforce" << i << "=" << c_force << "\n";
			}*/
			
		} //end fi
		if (pointradius(i) != 0)
		{
			multimastercontact(i)->SetLen(0);
			for (int fi=1; fi <= tempmasternum.Length(); fi++)
			{
				if (tempmasternum(fi) != 0) multimastercontact(i)->Add(tempmasternum(fi));
			}

		}
	}

	//if (nlerror > 0.01) UO() << "nlerror=" << nlerror << "\n";

	nlstepcnt++;
	TMStopTimer(25);

	return nlerror;
}

//vor SPH
//double GeneralContact2D::PostNewtonStep(double t) 
//{
//	TMStartTimer(25);
//
//	if (contactmode == 1 || contactmode == 2) //always contact
//	{
//		UO() << "ERROR: contact mode 1 and 2 not yet implemented for General Contact\n";
//		/*
//		rv = 0;
//		for (i = 1; i <= NSC(); i++)
//		{
//		if (!iscontact(i))
//		{
//		iscontact(i) = 1;
//		rv = 1;
//		}
//		}
//		return rv;
//		*/
//	}
//
//	int maxnlit = 3;
//	if (nlstepcnt > maxnlit) 
//	{
//		TMStopTimer(25);
//		return 0;
//	}
//	double nlerror = 0;
//
//	for (int i = 1; i <= NMC(); i++)
//	{
//		mastertoslave[i-1].SetLen(0);
//	}
//
//	for (int i = 1; i <= NSC(); i++)
//	{
//		Vector2D p, pp;
//		double gap;
//		int mnum, mlocind;
//
//		p = GetSlaveNodePos(i);
//
//		int n=1;
//		if (pointradius(i) == 0)
//		{
//			// get nearest segment or circle
//			GetNearestMasterSegment(i, p, mnum, mlocind, gap, pp); //gap < 0 means penetration
//			if (mnum == 0) n = 0;
//		}
//		else
//		{
//			GetNearestMasterSegments(i, p, pointradius(i), tempmasternum, templocind, tempdist, temp_pp);
//			n=tempmasternum.Length();
//		}
//
//		//if (tempmasternum.Length() && t > 2.83)
//		//	UO() << "master=" << tempmasternum << ", gaps=" << tempdist << "\n";
//
//
//		for (int fi=1; fi <= n; fi++)
//		{
//		//if (mnum != 0)
//		//{
//			if (pointradius(i) != 0)
//			{
//				mnum = tempmasternum(fi);
//				mlocind = templocind(fi);
//				gap = tempdist(fi);
//				pp = temp_pp(fi);
//			}
//
//			//UO() << "iscontact(" << i << ")=" << mnum << "\n";
//
//			masterlocind(mnum) = mlocind;
//
//			Vector2D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
//			Vector2D mpmid = GeomElem(mnum)->GetPos2D(ploc);
//			Vector2D n;
//			if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
//			else {n = p - mpmid; n.Normalize();}
//			Vector2D tt(-n.Y(), n.X());
//
//			double gapp;
//			Vector2D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
//			gapp = -((GetSlaveNodeVel(i) - mv) * n);      //project velocity to contact surface normal
//			double tvel = ((GetSlaveNodeVel(i) - mv) * tt); //project velocity tangential to surface
//
//			gap = (p-mpmid)*n - pointradius(i);
//
//
//			double c_force;
//			if (islagrange) 
//				c_force = XG(i);
//			else
//				c_force = GetContactForce(i, gap, gapp);
//
//			if (pointradius(i) != 0)
//			{
//				if (Find(mnum,*multimastercontact(i)))
//				{
//					iscontact(i) = mnum;
//				}
//				else
//				{
//					iscontact(i) = 0;
//					nlerror += 1;
//				}
//			}
//			else
//				if (iscontact(i) && mnum != iscontact(i)) 
//				{ 
//					iscontact(i) = mnum;
//					nlerror += 1;
//				}
//
//			if (gap < 0 && iscontact(i) == 0 && nlstepcnt < maxnlit)
//			{
//				iscontact(i) = mnum; //master element number
//				nlerror += fabs(gap);
//			}
//			else 	if (gap > master_adhesion_coefficient(mnum) && iscontact(i) != 0)
//			{
//				iscontact(i) = 0;
//				nlerror += fabs(c_force)/c_contact(i);
//			}
//
//			if (gapp_init(i) < fabs(gapp) && iscontact(i) != 0)
//			{
//				nlerror += fabs(gapp_init(i) - fabs(gapp))/Maximum(gapp_init(i),fabs(gapp));
//				gapp_init(i) = fabs(gapp);
//			}
//
//			if (isfriction && pointradius(i) == 0)
//			{
//				if (isstick(i))
//				{
//					Vector2D t_force;
//					Vector2D lastpp = GeomElem(isstick(i))->GetPos2D(laststickp(i));
//					GetTangentStickForce(i, tt, tvel, pp, lastpp, t_force);
//
//					//if (t >= 0.16 && t <= 0.1601) UO() << "slave" << i << ": gap=" << gap << ", cforce=" << c_force << ", tforce=" << t_force << "\n";
//
//					
//					if (t_force.Norm() > frictioncoeff_stick*c_force)
//					{
//						isstick(i) = 0;
//						nlerror += fabs(t_force.Norm() - frictioncoeff_stick*c_force)/c_contact(i);
//						if (abs(slipdir(i)) == 2) slipdir(i) = -slipdir(i) / 2;
//					}
//				}
//				else
//				{
//					//isstick(i) = mnum;
//					
//					if ((double)slipdir(i)*tvel < 0)
//					{
//						//assumed slip direction not equal tangential slip velocity -> try stick or other velocity:
//						isstick(i) = mnum;
//						nlerror += fabs(tvel);
//						slipdir(i) *= 2;
//					}
//				}
//			}
//
//			if (iscontact(i) != 0)
//			{
//				mastertoslave[mnum-1].Add(i);
//			}
//
//			if (pointradius(i) != 0)
//			{
//				tempmasternum(fi) = iscontact(i);
//			}
//
//			/*
//			if (gap < 0.0 && (t >= 0.0))
//			{
//				UO() << "t" << t << ", gap" << i << "=" << gap << ", cforce=" << c_force << ", mnum=" << mnum;
//				UO() << ", iscontact=" << iscontact(i) << ", err=" << nlerror << "\n";
//				//UO() << "cforce" << i << "=" << c_force << "\n";
//			}*/
//			
//		} //end fi
//		if (pointradius(i) != 0)
//		{
//			multimastercontact(i)->SetLen(0);
//			for (int fi=1; fi <= tempmasternum.Length(); fi++)
//			{
//				if (tempmasternum(fi) != 0) multimastercontact(i)->Add(tempmasternum(fi));
//			}
//
//		}
//	}
//
//	//if (nlerror > 0.01) UO() << "nlerror=" << nlerror << "\n";
//
//	nlstepcnt++;
//	TMStopTimer(25);
//
//	return nlerror;
//}



void GeneralContact2D::EvalG(Vector& f, double t) 
{
	/*
	UO() << "cmode=" << contactmode << "\n";
	UO() << "islagrange=" << islagrange << "\n";
	UO() << "isfriction=" << isfriction << "\n";*/

	if (!islagrange) return;
	
	if (isfriction) UO() << "ERROR: GeneralContact2D:EvalG, lagrange contact for friction not implemented!\n";

	if (slaveelements.Length() == 0)
	{
		UO() << "ERROR: GeneralContact2D::EvalG, number of elements == 0\n"; return;
	}

	if (MaxIndex()<=3)
	{
		Vector2D p1, p2, p, pp;
		for (int i = 1; i <= NSC(); i++)
		{
			if (pointradius(i) != 0) UO() << "ERROR: GeneralContact2D:EvalG, Lagrange contact for circle-elements not implemented!\n";

			//UO() << "iscontact(" << i << ")=" << iscontact(i) << "\n";
			if (iscontact(i) == 0)
			{
				f(i) = XG(i);
			}
			else 
			{
				double gap;
				p = GetSlaveNodePos(i);

				//mastersegment number is iscontact(i)
				int mnum = iscontact(i);


				gap = GeomElem(mnum)->GetNearestPoint(p, masterlocind(mnum), pp) - pointradius(i);

				Vector2D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
				Vector2D mpmid = GeomElem(mnum)->GetPos2D(ploc);
				Vector2D n;
				if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
				else {n = p - mpmid; n.Normalize();}
				Vector2D tt(-n.Y(), n.X());

				double gapp;
				Vector2D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
				gapp = -((GetSlaveNodeVel(i) - mv) * n);      //project velocity to contact surface normal
				//double tvel = ((GetSlaveNodeVel(i) - mv) * tt); //project velocity tangential to surface

				gap = (p-mpmid)*n - pointradius(i);


				switch (contactmode)
				{
				case 0: //Hertian contact
					f(i) = GetContactForce(i, gap, gapp) - XG(i); //XG(i) is compression force, Hertzian contact model
					break;
				case 1: //tied, linear elastic
					f(i) = GetContactForce(i, gap, gapp) - XG(i);
					break;
				case 2: //tied, rigid, Lagrange multiplier
					f(i) = gapp;
					break;
				case 3: //linear contact model
					f(i) = GetContactForce(i, gap, gapp) - XG(i);
					//f(i) = -gap - XG(i)/c_contact(i); //better condition number???
					break;
				case 4: //Lagrange multiplier contact, velocity based
					f(i) = gapp;
					break;
				default: UO() << "ERROR: contactmode\n"; f(i) = 0;
				}
			}
		}
	}

};

int cqcnt = 0;
int cqcnt2 = 0;
int cqcnt3 = 0;
void GeneralContact2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc

	//if (cqcnt%1000000 == 0) UO() << "AddCq=" << cqcnt << ", cq2=" << cqcnt2 << ", cq3=" << cqcnt3 << "\n";
	//cqcnt++;

	//if (cqcnt==1) UO() << "Ellen=" << elements.Length() << "\n";

	if (locelemind <= NSC())
	{
		//SLAVE
		int i = locelemind; //slave elements

		if ((iscontact(i) && pointradius(i) == 0) || (pointradius(i) != 0 && multimastercontact(i)->Length() != 0))
		{
			//cqcnt2++;

			double sign = 1;

			Vector2D p, pp;
			double gap;
			p = GetSlaveNodePos(i);

			int n = 1;
			if (pointradius(i) != 0) n = multimastercontact(i)->Length();
			for (int mi = 1; mi <= n; mi++)
			{

				//mastersegment number is iscontact(i)
				int mnum;
				if (pointradius(i) == 0) mnum = iscontact(i);
				else mnum = multimastercontact(i)->Get(mi);

				gap = GeomElem(mnum)->GetNearestPoint(p, masterlocind(mnum), pp) - pointradius(i);

				Vector2D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
				Vector2D mpmid = GeomElem(mnum)->GetPos2D(ploc);
				Vector2D n;
				if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
				else {n = p - mpmid; n.Normalize();}
				Vector2D tt(-n.Y(), n.X());

				double gapp;
				//Vector2D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
				Vector2D mv = GeomElem(mnum)->GetVel2D(ploc); //alternative, not so good

				Vector2D relv(GetSlaveNodeVel(i) - mv);
				gapp = -(relv * n);      //project velocity to contact surface normal
				double tvel = (relv * tt); //project velocity tangential to surface

				gap = (p-mpmid)*n - pointradius(i);


				double c_force;
				if (islagrange) 
					c_force = XG(i);
				else
					c_force = GetContactForce(i, gap, gapp);

				Vector2D totalforce(sign*n.X()*c_force, sign*n.Y()*c_force);

				if (isfriction)
				{
					if (isstick(i))
					{
						Vector2D t_force;
						Vector2D lastpp = GeomElem(isstick(i))->GetPos2D(laststickp(i));
						GetTangentStickForce(i, tt, tvel, pp, lastpp, t_force);

						totalforce += sign*t_force;

					}
					else
					{
						totalforce += -sign*slipdir(i)*frictioncoeff*fabs(c_force)*tt;
					}
				}

				if (slaveNODEmode)
				{
					int locn = locnodenums(i);
					double fact = 1.5;
					if (locn >= 1 && locn <= 4) fact = 3./4.;
					fact  = 1;
					totalforce *= fact; //correction of quadratic elements, //sign of totalforce: -f is incorporated in AddNodedPosdqTLambda

					GetSlaveBody(i).AddNodedPosdqTLambda(locnodenums(i),totalforce,f); 
				}
				else
				{
					GetSlaveBody(i).GetdPosdqT(loccoords(i),dpdq);
					//UO() << "dpdq" << i << "=" << dpdq << "\n";
					//UO() << "force=" << totalforce << "\n";


					for (int j=1; j <= f.Length(); j++)
					{
						f(j) -= -(dpdq(j,1)*totalforce.X()+dpdq(j,2)*totalforce.Y());
					}
				}
			}
		}

	}
	else
	{
		//MASTER:
		int j = locelemind - NSC(); //master elements

		if (GetMasterElnum(j))
		{
			double sign = -1;

			Vector2D p, pp;
			double gap;
			for (int k = 1; k <= mastertoslave[j-1].Length(); k++)
			{
				//cqcnt3++;
				int i = mastertoslave[j-1](k);

				p = GetSlaveNodePos(i);

				//mastersegment number is j
				int mnum = j;
				//mnum = iscontact(i);

				gap = GeomElem(mnum)->GetNearestPoint(p, masterlocind(mnum), pp) - pointradius(i);

				Vector2D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
				Vector2D mpmid = GeomElem(mnum)->GetPos2D(ploc);
				Vector2D n;
				if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
				else {n = p - mpmid; n.Normalize();}
				Vector2D tt(-n.Y(), n.X());

				double gapp;
				//Vector2D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
				Vector2D mv = GeomElem(mnum)->GetVel2D(ploc); //alternative, not so good
				Vector2D relv(GetSlaveNodeVel(i) - mv);
				gapp = -(relv * n);      //project velocity to contact surface normal
				double tvel = (relv * tt); //project velocity tangential to surface

				gap = (p-mpmid)*n - pointradius(i);


				/*
				if (gap < 1e-6 && GetMBS()->GetTime() > 0.001)
				{
				UO() << "master" << mnum << ": gap=" << gap << ", ploc=" << ploc << ", oldgap=" << oldgap << "\n";
				}*/

				double c_force;
				if (islagrange) 
					c_force = XG(i);
				else
					c_force = GetContactForce(i, gap, gapp);

				if (IsNaN(c_force))
				{
					UO() << "ERROR: c_force is NAN\n";
					UO() << "mnum=" << mnum << ", p=" << p << ", pp=" << pp << ", mpmid=" << mpmid << "\n   n=" << n << ", pointrad=" << pointradius(i) << "\n";
				}
				if (IsNaN(gap))
				{
					UO() << "ERROR: gap is NAN\n";
					UO() << "mnum=" << mnum << ", p=" << p << ", pp=" << pp << ", mpmid=" << mpmid << "\n   n=" << n << ", pointrad=" << pointradius(i) << "\n";
				}

				Vector2D totalforce(sign*n.X()*c_force, sign*n.Y()*c_force);
				//totalforce = Vector2D(0.,0.);

				if (isfriction)
				{
					//UO() << "ERROR:friction\n";
					if (isstick(i))
					{
						Vector2D t_force;
						Vector2D lastpp = GeomElem(isstick(i))->GetPos2D(laststickp(i));
						GetTangentStickForce(i, tt, tvel, pp, lastpp, t_force);

						totalforce += sign*t_force;

					}
					else
					{
						totalforce += -sign*slipdir(i)*frictioncoeff*fabs(c_force)*tt;
					}
				}

				GetMasterBody(j).GetdPosdqT(ploc,dpdq);

				for (int k=1; k <= f.Length(); k++)
				{
					f(k) -= -(dpdq(k,1)*totalforce.X()+dpdq(k,2)*totalforce.Y());
				}
			}
		}

	}
};

Vector3D GeneralContact2D::GetRefPosD() const 
{
	Vector2D v;
	if (NSC() != 0)
	{
		v = GetSlaveNodePosD(1);
	}
	return Vector3D(v.X(),v.Y(),0);
}

void GeneralContact2D::DrawElement() 
{
	Vector2D p1, p2, p, pp;
	for (int i = 1; i <= NSC(); i++)
	{
		if (GetMBS()->GetIOption(112))
		{
			p = GetSlaveNodePosD(i);

			//draw contact node:
			mbs->SetColor(GetCol());
			mbs->DrawSphere(GetSlaveBody(i).ToP3D(p),0.5*draw_dim.X(),4);

			/*
			//mastersegment number is iscontact(i)
			int mnum = iscontact(i);
			if (mnum && 0)
			{
				GeomElem(mnum)->GetNearestPoint(p, masterlocind(mnum), pp);

				mbs->DrawSphere(GetSlaveBody(i).ToP3D(pp),1.*draw_dim.X(),3);
			}*/
		}
	}
	if (GetMBS()->GetIOption(112))
	{
		for (int i = 1; i <= NMC(); i++)
		{
			//if (GetMBS()->GetIOption(112))
			if (masterisactive(i) == 0 || GetMasterIsActiveCondition(GetMBS()->GetDrawTime(), i))
			{
				GeomElem(i)->DrawYourself();
			}
		}
	}
};

	// AP 2013-01-11: removed source code to cpp file, otherwise FinishContactDefinition cannot be called in model file
void GeneralContact2D::FinishContactDefinition()
{
	x_init = Vector(IS()); //number of normal contact lagrange multipliers == number of slave nodes

	int nnc = NSC(); //number of normal contacts
	iscontact.SetLen(nnc);
	gapp_init.SetLen(nnc);
	gapp_init_min.SetLen(nnc);

	isstick.SetLen(nnc);
	islaststick.SetLen(nnc);
	laststickp.SetLen(nnc);
	slipdir.SetLen(nnc);
	multimastercontact.SetLen(nnc); //only elements set where pointradius != 0

	for (int i=1; i <= NSC(); i++)
	{
		iscontact(i) = 0;
		isstick(i) = 0;
		islaststick(i) = 0;
		laststickp(i) = Vector2D(0.,0.);
		slipdir(i) = 1;

		gapp_init_min(i) = 1e-6; //1e-2 for hourglass
		gapp_init(i) = gapp_init_min(i);

		if (pointradius(i) != 0)
		{
			multimastercontact(i) = new IVector(2);
		}
	}
	elements.SetLen(0);
	for (int i = 1; i <= NSC(); i++)
	{	
		elements.Add(slaveelements(i));
	}
	for (int i = 1; i <= NMC(); i++)
	{	
		elements.Add(masterelements(i)->GetElnum());
	}

	//for element jacobian:
	elements_nodouble.SetLen(0);
	for (int i = 1; i <= NSC(); i++)
	{	
		if (!Find(slaveelements(i), elements_nodouble))
		{
			elements_nodouble.Add(slaveelements(i));
		}
	}

	for (int i = 1; i <= NMC(); i++)
	{	
		if (!Find(masterelements(i)->GetElnum(), elements_nodouble))
		{
			elements_nodouble.Add(masterelements(i)->GetElnum());
		}
	}

	mastertoslave = new IVector[NMC()]();
	for (int i = 0; i < NMC(); i++)
	{
		mastertoslave[i].Init();
	}

	tempitemlist.SetLen(NMC());
	for (int i = 1; i <= NMC(); i++)
		tempitemlist(i) = 0;
}