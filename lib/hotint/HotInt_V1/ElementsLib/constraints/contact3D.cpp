//#**************************************************************
//#
//# filename:             contact3D.cpp
//#
//# author:               Gerstmayr Johannes, Peter Gruber
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
#include "body3d.h"
#include "geomelements.h"
#include "contact3d.h"
#include "femathhelperfunctions.h"
#include "contact_globalvariables.h"


int contact_cnt_gap = 0;
int contact_cnt_mastersearch = 0;
int contact_cnt_itemlist = 0;
int contact_cnt_foundmaster = 0;
int contact_cnt_list2 = 0;
int contact_cnt_list3 = 0;
int contact_cnt_list4 = 0;

int testcnt = 1000;
int last_buildsearchtree = 1000;

const int useconstantmastersearchtree = 0;
const int discstepcorrection = 1; //correct energy lost during discontinuous step (contact switching)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GeneralContact3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GeneralContact3D::BuildSearchTrees()
{

	if (!mastersearchtreeinitialized)
	{
		//+++++++++++++++++++++++++++++
		//master segments:
		Box3D box;
		box.Clear();
		for (int i=1; i <= NMC(); i++)
		{
			Box3D pbox;
			GeomElement* ge = GeomElem(i);
			pbox = ge->GetBoundingBox();
			box.Add(pbox);
		}
		//UO() << "STbox=" << box.PMin() << "-" << box.PMax() << "\n";
		box.Increase(bordersize);
		mastersearchtreebox = box; //this box is currently not rebuilt-->if objects move outside initial box, search will be slower!

		if (useconstantmastersearchtree)
		{
			//non-moving-objects:
			constantmastersearchtree.ResetSearchTree((int)(1.5*searchtreeix), (int)(1.5*searchtreeiy), (int)(1.5*searchtreeiz), mastersearchtreebox); //define number of cubes and maximum size of search tree
			for (int i=1; i <= NMC(); i++)
			{
				if (GeomElem(i)->GetElnum() == 0)
				{
					Box3D pbox = GeomElem(i)->GetBoundingBox();

					pbox.Increase(bordersize);
					constantmastersearchtree.AddItem(pbox,i);
				}
			}
		}
		mastersearchtreeinitialized = 1;
	}

	if (last_buildsearchtree++ < 0) return; //save time when not updating search trees so often

	last_buildsearchtree = 0;

	TMStartTimer(24); //build search trees

	mastersearchtree.ResetSearchTree(searchtreeix, searchtreeiy, searchtreeiz, mastersearchtreebox); //define number of cubes and maximum size of search tree


	for (int i=1; i <= NMC(); i++)
	{
		if (GeomElem(i)->GetElnum() != 0 || !useconstantmastersearchtree)
		{
			Box3D pbox = GeomElem(i)->GetBoundingBox();

			pbox.Increase(bordersize);
			mastersearchtree.AddItem(pbox,i);
		}
	}
	TMStopTimer(24); //build search trees
}

const Body3D& GeneralContact3D::GetSlaveBody(int i) const 
{
	return (const Body3D&)mbs->GetElement(slaveelements(i));
}

Body3D& GeneralContact3D::GetSlaveBody(int i) 
{
	return (Body3D&)mbs->GetElement(slaveelements(i));
}

//$ YV 2012-12-13: replaced GeomElem(i)->GetBody3D() by GeomElem(i)->GetElement() here and below
const Body3D& GeneralContact3D::GetMasterBody(int i) const
{
	return (Body3D&)GeomElem(i)->GetElement();
}

Body3D& GeneralContact3D::GetMasterBody(int i) 
{
	return (Body3D&)GeomElem(i)->GetElement();
}

Vector3D GeneralContact3D::GetSlaveNodePos(int i) const 
{
	if (slaveNODEmode == 0)
	{
		return GetSlaveBody(i).GetPos(loccoords(i));
	}
	else
	{
		return GetSlaveBody(i).GetNodePos(locnodenums(i));
	}
}

Vector3D GeneralContact3D::GetSlaveNodePosD(int i) const 
{
	if (slaveNODEmode == 0)
	{
		return GetSlaveBody(i).GetPosD(loccoords(i));
	}
	else
	{
		return GetSlaveBody(i).GetNodePosD(locnodenums(i));
	}
}

Vector3D GeneralContact3D::GetSlaveNodeVel(int i) const 
{
	if (slaveNODEmode == 0)
	{
		return GetSlaveBody(i).GetVel(loccoords(i));
	}
	else
	{
		return GetSlaveBody(i).GetNodeVel(locnodenums(i));
	}
}

Vector3D GeneralContact3D::GetMasterLocPos(const Vector3D& pp, int j, int locind) const
{
	return GeomElem(j)->GetLocPos(locind, pp);
}

//approximate velocity at segment point pglob due linear interpolation
Vector3D GeneralContact3D::GetMasterSegVel(const Vector3D& pglob, int j, int locind) const 
{
	return GeomElem(j)->GetVel(GeomElem(j)->GetLocPos(locind, pglob));
}

//get mindist and projected point pp; only for points, spheres return several elements!!!
void GeneralContact3D::GetNearestMasterSegment(int slavenum, const Vector3D& p, int& masternum, 
																							 int& locind, double& dist, Vector3D& pp)
{
	Box3D box;
	box.Add(p);
	//box.Increase(pointradius(slavenum));

	mastersearchtree.GetItemsInBox(box, itemlist);
	if (useconstantmastersearchtree) constantmastersearchtree.AddItemsInBox(box, itemlist);

	Vector3D pf;

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

	for (int i = 1; i <= itemlist2.Length(); i++)
	{
		int j = itemlist2(i);
		tempitemlist(j) = 0;

		if (masterbodyindex(j) < slavebodyindex(slavenum) &&
				(masterisactive(j) == 0 || mathfunctions(masterisactive(j))->CheckCondition(GetMBS()->GetTime())) )
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

//circle radius is radius of sphere!!!
void GeneralContact3D::GetNearestMasterSegments(int slavenum, const Vector3D& p, double circlerad, TArray<int>& a_masternum, 
		TArray<int>& a_locind, TArray<double>& a_dist, TArray<Vector3D>& a_pp)
{
	contact_cnt_mastersearch++; //400.000/sec //50%

	//TMStartTimer(28); //60% of GetNearestMasterSegs
	Box3D box;
	box.Add(p);
	box.Increase(circlerad);

	//int ind[6]; //indices for box
	//mastersearchtree.GetBoxIndizes(box, ind); //8%
	//mastersearchtree.GetItemsInBox(ind, itemlist); //10.8%

	mastersearchtree.GetItemsInBox(box, itemlist); //2.500.000 Elements/sec
	if (useconstantmastersearchtree) constantmastersearchtree.AddItemsInBox(box, itemlist);

	Vector3D pf;
	double d;
	int flocind = 0;
	double gap;

	int	locind = 0; 
	int masternum = 0;
	double dist = 1e20;
	Vector3D pp;

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
	contact_cnt_itemlist += itemlist2.Length();
//itemlist = 33, itemlist2 = 10.027, list2=4.75744, list3=0.213187, foundmaster = 0.0111501
//n_gap=7.34539e+006, n_ms=3.9368e+007, itemlist=41.9261, itemlist2=13.5976, foundmaster=0.071787, list2=6.87892, list3=1.0026, 

	//TMStopTimer(28); //59% time spent up to here
	Vector3D pc;
	double rz;

	for (int i = 1; i <= itemlist2.Length(); i++)
	{
		int j = itemlist2(i);
		tempitemlist(j) = 0;

		if (masterbodyindex(j) < slavebodyindex(slavenum) &&
				(masterisactive(j) == 0 || mathfunctions(masterisactive(j))->CheckCondition(GetMBS()->GetTime())) )
		{
			//contact_cnt_list2++; //47% of items get here
			if (GeomElem(j)->GetType() == TGeomSphere3D && circlerad != 0)
			{
			//TMStartTimer(29);
				//contact_cnt_list3++;
				locind  = 1; //only one index possible in Sphere

				pc = GeomElem(j)->GetElement().GetRefPos();		//center of sphere, global; ACCELERATE: getrefpos directly?
				rz = ((GeomSphere3D*)GeomElem(j))->GetRadius();
				pf = p-pc;										//vector from center pc of sphere to p
				gap = pf.Norm2();							//distance of pc and p

				//gap = sqrt(gap);
				//d = gap-rz-circlerad;
				//if ((d < 0) != (gap < Sqr(rz+circlerad))) UO() << "problem in contact\n";

				//if (d < 0)  
				if (gap < Sqr(rz+circlerad))
				{
					//contact_cnt_list4++;
					gap = sqrt(gap);
					if (gap != 0) pf *= (rz/gap);	//pf=relative vector from center to surface of sphere
					pf += pc;											//pf=global projected point at surface of sphere
					d = gap-rz-circlerad;
					gap = d;											//distance between p and projected point pf
				}
				else
				{
					d = 1; //so that d > 0 !!!
				}
			//TMStopTimer(29);//30% of total time here
			}
			else
			{
				gap = GeomElem(j)->GetNearestPoint(p, locind, pf) - circlerad; //gap value??
				d = Dist(p,pf) - circlerad;

				if (fabs(d) < fabs(dist) && fabs(d) < contactmaxdist) //only for case when circlerad==0
				{
					dist = d;
					masternum = j;
					pp = pf;
					flocind = locind;
				}
			}

			if (circlerad != 0 && d < 0	&& fabs(d) < contactmaxdist)  
			{
				//0.1% of items get here
				//this is practically done for every master element
//				if (!GetMBS()->IsJacobianComputation() && 0)
//					UO() << "GeomElem" << j << ", gap=" << d << "\n";

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
	contact_cnt_foundmaster += a_masternum.Length();
	//otherwise leave arrays empty
}

void GeneralContact3D::GetTangentStickForce(int i, const Vector3D& t, double tvel, const Vector3D& pp, const Vector3D& lastpp, Vector3D& tforce)
{
  double friction_penalty = 0.1*c_contact(i);

	tforce = 0;
	if (islaststick(i) && isstick(i)) tforce = friction_penalty*(lastpp-pp);
	//else tforce = (-friction_penalty*tvel)*t;
}

double GeneralContact3D::ConvertRestitution(double restcoeff)
{
	double x = restcoeff;
	return 0.02842673968087363+1.235462852658231*pow(x,7)-27.76611509389485*pow(x,4)-8.758510573837406*pow(x,6)+22.67054845547086*pow(x,5)+16.29462160441537*pow(x,3)
		-4.061556592190019*pow(x,2)+1.357122607696947*x;
}

double GeneralContact3D::GetContactForce(int i, double gap, double gapp)
{
	//double e = 0.95; //coeff of restitution; 0.95 ususally, with ambrosio
	double e = restitution_coeff; //coeff of restitution; 0.95 ususally
	double nh = hertzian_contact_param; //coeff for Hertzian contact; 1 works good

	double gi = fabs(gapp_init(i));
	if (fabs(gapp) > fabs(gapp_init(i))) gi = fabs(gapp);

	switch (contactmode)
	{
	case 0: //Hertian contact
		return  -Sgn(gap)*c_contact(i)*(1.+0.75*(1.-Sqr(e))*(-gapp)/gi)*pow(fabs(gap),nh); //original, XG(i) is compression force, Hertzian contact model
		//following model has been used in some numerical examples:
		//return  -Sgn(gap)*c_contact(i)*(1.+8*(1.-Sqr(e))*gapp/gi)*pow(fabs(gap),nh); //XG(i) is compression force, Hertzian contact model
		break;
	case 1: //tied, linear elastic
		return  -c_contact(i)*gap;
		break;
	case 3: //linear contact model
		return  -c_contact(i)*gap;
		break;
	case 5: //Model which gives correct restitution coefficient from 0..1
		return  -Sgn(gap)*c_contact(i)*(1.+0.75*(1.-Sqr(e))/(1.-Sqr(1.-e))*(-gapp)/gi)*pow(fabs(gap),nh); //XG(i) is compression force, Hertzian contact model
		break;
	default: UO() << "ERROR: contactmode in GetContactForce\n"; return 0;
	}
}

void GeneralContact3D::PostprocessingStep() 
{
	//TMStartTimer(24);
	for (int i = 1; i <= NMC(); i++)
	{
		mastertoslave_last[i-1] = mastertoslave[i-1];
		mastertoslave_actlast[i-1] = mastertoslave[i-1];
	}

	for (int i = 1; i <= NSC(); i++)
	{
		if (pointradius(i) != 0) 
		{
			*multimastercontact_last(i) = *multimastercontact(i);
			*multimastercontact_actlast(i) = *multimastercontact(i);
			if (multimastercontact_last(i)->Length() == 0) gapp_init(i) = gapp_init_min(i); //reset should be done earlier for each contact pair!!!
		}
	}

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

						Vector3D pglob = GetSlaveNodePos(i); 

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

						Vector3D pglob = GetSlaveNodePos(i);

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

//ofstream gapout("..\\..\\output\\gap.txt"); //Ö

double GeneralContact3D::PostNewtonStep(double t) 
{
	TMStartTimer(25);

	if (contactmode == 1 || contactmode == 2) //always contact
	{
		UO() << "ERROR: contact mode 1 and 2 not yet implemented for General Contact\n";
	}

	int maxnlit = 4; //2 works fine, but ignores many special cases
	if (nlstepcnt > maxnlit) 
	{
		TMStopTimer(25);
		return 0;
	}
	double nlerror = 0;

	for (int i = 1; i <= NMC(); i++)
	{
		mastertoslave[i-1].SetLen(0);
	}

	for (int i = 1; i <= NSC(); i++) //i is index of slave element
	{
		Vector3D p, pp;
		double gap;
		int mnum, mlocind;

		p = GetSlaveNodePos(i);

		int nseg=1;
		if (pointradius(i) == 0)
		{
			GetNearestMasterSegment(i, p, mnum, mlocind, gap, pp); //gap < 0 means penetration; mnum is master element
			if (mnum == 0) nseg = 0;
		}
		else
		{
			TMStartTimer(27);
			GetNearestMasterSegments(i, p, pointradius(i), tempmasternum, templocind, tempdist, temp_pp);
			nseg=tempmasternum.Length();
			TMStopTimer(27);
		}

		//?????????????
		//Add multimastercontact_actlast elements somehow to tempmasternum_actlast list ...
		//then switching back from contact will be implemented correctly

		for (int fi=1; fi <= nseg; fi++)
		{
		//if (mnum != 0)
		//{
			if (pointradius(i) != 0)
			{
				mnum = tempmasternum(fi); //mnum is master element
				mlocind = templocind(fi);
				gap = tempdist(fi);
				pp = temp_pp(fi);
			}

			double gapp, tvel, forcefact, forceadd;
			Vector3D t1, t2, n, ploc;

			//gapp is d/dt(gap)
			ComputeGap(i, mnum, gap, gapp, p, pp, ploc, t1, t2, n, tvel, forcefact, forceadd, 1); //called from PostNewtonStep

			double c_force;
			if (islagrange) 
				c_force = XG(i);
			else
				c_force = GetContactForce(i, gap, gapp);

			c_force *= forcefact;
			c_force += forceadd;

			if (pointradius(i) != 0)
			{
				if (Find(mnum,*multimastercontact(i))) //compare with last nonlinear iteration / assumption of contact!
				{
					iscontact(i) = mnum; //this variable is reused in multimastercontact mode ....!
					//if (GetMBS()->GetStepSize() > 1e-7) nlerror = -1;//GetMBS()->GetStepSize(); //adaptive stepsize for contact
				}
				else
				{
					iscontact(i) = 0;
					nlerror += 1;
				}
			}
			else if (iscontact(i) && mnum != iscontact(i)) 
			{ 
				iscontact(i) = mnum;
				nlerror += 1;
			}

			if (gap < 0 && iscontact(i) == 0 && nlstepcnt < maxnlit)
			{
				iscontact(i) = mnum; //master element number
				nlerror += fabs(gap);
			}
			else 	if (gap > 0 && iscontact(i) != 0)
			{
				iscontact(i) = 0;
				nlerror += fabs(c_force)/c_contact(i);
			}

			if (fabs(gapp_init(i)) < fabs(gapp) && iscontact(i) != 0)
			{
				nlerror += fabs(fabs(gapp_init(i)) - fabs(gapp))/Maximum(fabs(gapp_init(i)),fabs(gapp));
				gapp_init(i) = fabs(gapp);
			}

			if (isfriction && pointradius(i) == 0)
			{
				if (isstick(i))
				{
					UO() << "Contact3D::not implemented!!!\n";
					Vector3D t_force;
					Vector3D lastpp = GeomElem(isstick(i))->GetPos(laststickp(i));
					GetTangentStickForce(i, t1, tvel, pp, lastpp, t_force);

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
			
		} //end fi
		if (pointradius(i) != 0)
		{
			//check if contact was active in last step, but is not active in current step:
			for (int fi=1; fi <= multimastercontact(i)->Length(); fi++)
			{
				if (!Find(multimastercontact(i)->Get(fi),tempmasternum))
					nlerror += 1;
			}

			multimastercontact(i)->SetLen(0);
			for (int fi=1; fi <= tempmasternum.Length(); fi++)
			{
				if (tempmasternum(fi) != 0) multimastercontact(i)->Add(tempmasternum(fi));
			}

			//multimastercontact_actlast should contain all actual contacts and contacts which were active in last time step
			*multimastercontact_actlast(i) = *multimastercontact(i);
			for (int fi = 1; fi <= multimastercontact_last(i)->Length(); fi++)
			{
				if (!Find(multimastercontact_last(i)->Get(fi),*multimastercontact_actlast(i))) 
				{
					multimastercontact_actlast(i)->Add(multimastercontact_last(i)->Get(fi));

					//check gap in this elements as well:
					double gapp, tvel, forcefact, forceadd;
					Vector3D t1, t2, n, ploc;
					ComputeGap(i, multimastercontact_last(i)->Get(fi), gap, gapp, p, pp, ploc, t1, t2, n, tvel, forcefact, forceadd, 1); //called from PostNewtonStep
				}
			}
			//if (multimastercontact(i)->Length() || multimastercontact_last(i)->Length()) 
			//UO() << "MMC=" << *multimastercontact(i) << ", mmc_last=" << *multimastercontact_last(i) << ", mmc_actlast=" << *multimastercontact_actlast(i) << "\n";
		}
	}

	for (int mnum = 1; mnum <= NMC(); mnum++) //update master to slave actual list
	{
		mastertoslave_actlast[mnum-1] = mastertoslave[mnum-1];

		for (int fi = 1; fi <= mastertoslave_last[mnum-1].Length(); fi++)
		{
			if (!Find(mastertoslave_last[mnum-1].Get(fi),mastertoslave_actlast[mnum-1])) 
			{
				mastertoslave_actlast[mnum-1].Add(mastertoslave_last[mnum-1].Get(fi));
				//UO() << "contact M" << mnum << "-" << mastertoslave_last[mnum-1].Get(fi) << " active before\n";
			}
		}
	}

	//if (nlerror > 0.01) UO() << "nlerror=" << nlerror << "\n";

	nlstepcnt++;
	TMStopTimer(25);

	return nlerror;
}


void GeneralContact3D::EvalG(Vector& f, double t) 
{
	/*
	UO() << "cmode=" << contactmode << "\n";
	UO() << "islagrange=" << islagrange << "\n";
	UO() << "isfriction=" << isfriction << "\n";*/

	if (!islagrange) return;
	
	if (isfriction) UO() << "ERROR: GeneralContact3D:EvalG, lagrange contact for friction not implemented!\n";

	if (slaveelements.Length() == 0)
	{
		UO() << "ERROR: GeneralContact3D::EvalG, number of elements == 0\n"; return;
	}

	if (MaxIndex()<=3)
	{
		Vector3D p, pp;
		for (int i = 1; i <= NSC(); i++)
		{
			if (pointradius(i) != 0) UO() << "ERROR: GeneralContact3D:EvalG, Lagrange contact for circle-elements not implemented!\n";

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

				Vector3D ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
				Vector3D mpmid = GeomElem(mnum)->GetPos(ploc);
				Vector3D n;
				if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
				else {n = p - mpmid; n.Normalize();}

				Vector3D t1 = GeomElem(mnum)->GetTangent(masterlocind(mnum), ploc);
				Vector3D t2 = n.Cross(t1);

				double gapp;
				Vector3D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
				gapp = ((GetSlaveNodeVel(i) - mv) * n);      //project velocity to contact surface normal
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
				case 5: //corrected restitution
					f(i) = GetContactForce(i, gap, gapp) - XG(i);
					break;
				default: UO() << "ERROR: contactmode\n"; f(i) = 0;
				}
			}
		}
	}

};

void GeneralContact3D::ComputeGap(int i, int mnum, double& gap, double& gapp, Vector3D& p, Vector3D& pp, 
																	Vector3D& ploc, Vector3D& t1, Vector3D& t2, Vector3D& n, double& tvel, double& forcefact, double& forceadd, int nlstep = 0)
{
contact_cnt_gap++;

TMStartTimer(26); //gap
	gap = GeomElem(mnum)->GetNearestPoint(p, masterlocind(mnum), pp) - pointradius(i);

	ploc = GetMasterLocPos(pp, mnum, masterlocind(mnum));
	Vector3D mpmid = GeomElem(mnum)->GetPos(ploc);

	if (pointradius(i) == 0) n = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc); //normal, outwards
	else {n = p - mpmid; n.Normalize();}
	t1 = GeomElem(mnum)->GetTangent(masterlocind(mnum), ploc);
	t2 = n.Cross(t1);

	Vector3D mv = GetMasterSegVel(pp, mnum, masterlocind(mnum)); //master segment approximate velocity
	gapp = ((GetSlaveNodeVel(i) - mv) * n);      //project velocity to contact surface normal

	//project velocity tangential to surface
	tvel = sqrt(Sqr((GetSlaveNodeVel(i) - mv) * t1) + Sqr((GetSlaveNodeVel(i) - mv) * t2)); 

	gap = (p-mpmid)*n - pointradius(i); //better approximation of gap


	//if (GetMBS()->GetTime() >= 0.2202 && nlstep) UO() << "t=" << GetMBS()->GetTime() << ", gap=" << gap << ", gapp=" << gapp << ", dt=" << GetMBS()->GetStepSize() << "\n";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int iscontactA = 1;			//contact at beginning of time step
	int iscontactB = 1;			//contact at end of time step

	forcefact = 1;		//correction of force, because contact force is not applied for the whole time step
	forceadd = 0;

	if (!Find(mnum, *multimastercontact_last(i)))
		//if (!Find(i, mastertoslave_last[mnum-1])) //check if contact was active at beginning of step!
	{
		iscontactA = 0;
	}
	if (!Find(mnum, *multimastercontact(i)))
		//if (!Find(i, mastertoslave[mnum-1])) //check if contact is active now!
	{
		iscontactB = 0;
	}


	int simplemode = 1;

	if (!simplemode)
	{
		//correction of force:
		if (iscontactA != iscontactB)
		{
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//compute gap0 and gapp0:
			double* oldptr = GetMBS()->GetXact().GetVecPtr();
			int oldlen = GetMBS()->GetXact().Length();
			GetMBS()->SetActState(GetMBS()->GetLastNLItSolVector());

			Vector3D p0 = GetSlaveNodePos(i);
			Vector3D pp0, ploc0, n0, t10, t20;

			double gap0 = GeomElem(mnum)->GetNearestPoint(p0, masterlocind(mnum), pp0) - pointradius(i);

			ploc0 = GetMasterLocPos(pp0, mnum, masterlocind(mnum));
			Vector3D mpmid0 = GeomElem(mnum)->GetPos(ploc0);

			if (pointradius(i) == 0) n0 = GeomElem(mnum)->GetNormal(masterlocind(mnum), ploc0); //normal, outwards
			else {n0 = p0 - mpmid0; n0.Normalize();}
			t10 = GeomElem(mnum)->GetTangent(masterlocind(mnum), ploc0);
			t20 = n0.Cross(t10);

			Vector3D mv0 = GetMasterSegVel(pp0, mnum, masterlocind(mnum)); //master segment approximate velocity
			double gapp0 = ((GetSlaveNodeVel(i) - mv0) * n0);      //project velocity to contact surface normal

			//project velocity tangential to surface
			double tvel0 = sqrt(Sqr((GetSlaveNodeVel(i) - mv0) * t10) + Sqr((GetSlaveNodeVel(i) - mv0) * t20)); 

			gap0 = (p0-mpmid0)*n0 - pointradius(i); //better approximation of gap

			GetMBS()->SetActState(oldptr, oldlen);
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			//UO() << "cA=" << iscontactA << ", cB=" << iscontactB << "\n";
			//gapp = -d/dt(gap)  .......
			double dt = GetMBS()->GetStepSize();
			double gaptot = fabs(gap0) + fabs(gap);
			if (dt != 0 && gaptot != 0)
			{
				double t2 = fabs(gap)/gaptot*dt; //estimated time since last switch
				double t1 = dt-t2;

				//t2 = gap/gapp;

				//UO() << "t2/dt=" << t2/dt << ",gap0=" << gap0 << ",gap=" << gap << ",gapp0=" << gapp0 << ",gapp=" << gapp << "\n";
				int conservingmode = 1;
				if (nlstep)
				{
					if (t1 < GetMBS()->GetStepRecommendation() && t1 < GetMBS()->GetStepSize() && fabs(gap) > 2.e-7) 
					{
						GetMBS()->SetStepRecommendation(t1);
						//if (t1 > 1e-6) GetMBS()->SetStepRecommendation(1e-6);

						//UO() << "contact rec dt=" << t1;// << "\n";	
						//UO() << "c" << iscontactA << iscontactB << ", t=" << GetMBS()->GetTime() << ", gap=" << gap << ", gapp=" << gapp << ", dt=" << GetMBS()->GetStepSize() << "\n";				
					}
				}

				if (!iscontactA && iscontactB)
				{
					if (conservingmode)
					{
						forcefact = 0;// 0
						//forceadd = -1.361207*c_contact(i)*(gap*(t2/dt)); //for first step right with forcefact = 0 ...

						//forceadd = GetContactForce(i, gap, gapp)*(t2/dt)*(1.+2.*Sqr(t1/dt)); //for first step right with forcefact = 0 ...
						forceadd = GetContactForce(i, gap, gapp)*(t2/dt); //for first step right with forcefact = 0 ...
					}
					else
					{
						forcefact = 1;// 0
						forceadd = 0;
					}
					//forceadd =  -c_contact(i)*gap0*(t1/dt);

					//forceadd  = -c_contact(i)*(-Sqr(t2)*dt*gapp0+Sqr(dt)*gap0+Cub(dt)*gapp0-Sqr(t2)*gap0)/Sqr(dt); //equation from acceleration
					//forceadd  = -c_contact(i)*t2*(dt*gap0+dt*t1*gapp0-t2*gap0)/Sqr(dt); //equation from velocities

					//if (!GetMBS()->IsJacobianComputation()) UO() << "t=" << GetMBS()->GetTime() << ", t2/dt=" << t2/dt << ", fact=" << forcefact << ", f=" << forceadd << ", gap=" << gap << ", gapp=" << gapp << "\n";
				}

				if (iscontactA && !iscontactB && gap != 0)
				{
					if (conservingmode)
					{
						forcefact = 0; //too big???
						//forceadd = -c_contact(i)*gap0*(t1/dt); //somehow right with forcefact = 0 ...
						//forceadd = c_contact(i)*gap*(t2/dt)*(1.+2.*Sqr(t1/dt)); //somehow right with forcefact = 0 ...
						forceadd = -0*GetContactForce(i, gap, gapp)*(t2/dt); //somehow right with forcefact = 0 ...

						//if (!GetMBS()->IsJacobianComputation()) UO() << "t=" << GetMBS()->GetTime() << ", t2/dt=" << t2/dt << ", fact=" << forcefact << ", f=" << forceadd << ", gap=" << gap << ", gapp=" << gapp << "\n";
					}
					else
					{
						forcefact = 0; //too big???
						forceadd = 0;
						gap = 0;
					}
				}
			} 
		}
	}
	else
	{
		forcefact = 1;
		forceadd = 0;
	}
TMStopTimer(26); //gap
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}


void GeneralContact3D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;

	Vector3D p, pp;
	double gap;

	if (locelemind <= NSC())
	{
		//SLAVE
		int i = locelemind; //slave elements

		//if (GetMBS()->GetTime() > 0.5704)			UO() << "t2=" << GetMBS()->GetTime() << "mmc=" << *multimastercontact(i) << ", mmc_last=" << *multimastercontact_last(i) << ", mmc_actlast=" << *multimastercontact_actlast(i) << "\n";
 
		if ((pointradius(i) == 0 && iscontact(i)) || (pointradius(i) != 0 && multimastercontact_actlast(i)->Length() != 0))
		{
			p = GetSlaveNodePos(i);

			int n = 1;
			if (pointradius(i) != 0) n = multimastercontact_actlast(i)->Length();
			for (int mi = 1; mi <= n; mi++)
			{

				//mastersegment number is iscontact(i)
				int mnum;
				if (pointradius(i) == 0) mnum = iscontact(i);
				else mnum = multimastercontact_actlast(i)->Get(mi);

				double gapp, tvel, forcefact, forceadd;
				Vector3D t1, t2, n, ploc;
 
				//gapp is d/dt(gap)
				ComputeGap(i, mnum, gap, gapp, p, pp, ploc, t1, t2, n, tvel, forcefact, forceadd);

				double c_force;
				if (islagrange) 
					c_force = XG(i);
				else
					c_force = GetContactForce(i, gap, gapp);

				c_force *= forcefact;
				c_force += forceadd;

				Vector3D totalforce(sign*n.X()*c_force, sign*n.Y()*c_force, sign*n.Z()*c_force);
				
				/*
				if (gap < 0.0 && (GetMBS()->GetTime() >= 0.0))
				{
					UO() << "gap" << i << "=" << gap << ", force=" << totalforce << "\n";
				}*/

				if (isfriction)
				{
					UO() << "Contact3D:: Tangent stick force must be corrected for friction!\n";
					if (isstick(i))
					{
						Vector3D t_force;
						Vector3D lastpp = GeomElem(isstick(i))->GetPos(laststickp(i));
						GetTangentStickForce(i, t1, tvel, pp, lastpp, t_force); //not correct!!!

						totalforce += sign*t_force;

					}
					else
					{
						//slipdir->slipvector2D !!!
						//totalforce += -sign*slipdir(i)*frictioncoeff*fabs(c_force)*(t1...t2); 
					}
				}

				if (slaveNODEmode)
				{
					int locn = locnodenums(i);
					double fact;
					/*
					fact = 1.5;
					if (locn >= 1 && locn <= 4) fact = 3./4.;
					*/
					fact  = 1;
					totalforce *= -fact; //correction of quadratic elements, //negative sign: needs to be subtracted from f!!!!

					GetSlaveBody(i).AddNodedPosdqTLambda(locnodenums(i),totalforce,f); 
				}
				else
				{
					GetSlaveBody(i).GetdPosdqT(loccoords(i),dpdq);
					//UO() << "dpdq" << i << "=" << dpdq << "\n";
					//UO() << "force=" << totalforce << "\n";
					//UO() << "loccoords=" << loccoords(i) << "\n";


					for (int j=1; j <= f.Length(); j++)
					{
						f(j) -= -(dpdq(j,1)*totalforce.X()+dpdq(j,2)*totalforce.Y()+dpdq(j,3)*totalforce.Z());
					}
				}
			}
		}
	}
	else
	{
		//MASTER:
		int j = locelemind - NSC(); //master elements
		sign = -1;

		if (GetMasterElnum(j))
		{
			for (int k = 1; k <= mastertoslave_actlast[j-1].Length(); k++)
			{
				int i = mastertoslave_actlast[j-1](k);

				p = GetSlaveNodePos(i);

				//mastersegment number is j
				int mnum = j;
				//mnum = iscontact(i);


				double gapp, tvel, forcefact, forceadd;
				Vector3D t1, t2, n, ploc;
 
				ComputeGap(i, mnum, gap, gapp, p, pp, ploc, t1, t2, n, tvel, forcefact, forceadd);

				double c_force;
				if (islagrange) 
					c_force = XG(i);
				else
					c_force = GetContactForce(i, gap, gapp);

				c_force *= forcefact;
				c_force += forceadd;

				if (IsNaN(c_force))
					UO() << "ERROR: c_force is NAN\n";
				if (IsNaN(gap))
				{
					UO() << "ERROR: gap is NAN\n";
					UO() << "mnum=" << mnum << ", p=" << p << ", pp=" << pp << "\n";
				}

				Vector3D totalforce(sign*n.X()*c_force, sign*n.Y()*c_force, sign*n.Z()*c_force);

				if (isfriction)
				{
					UO() << "Contact3D::friction, master\n";
					if (isstick(i))
					{
						Vector3D t_force;
						Vector3D lastpp = GeomElem(isstick(i))->GetPos(laststickp(i));
						GetTangentStickForce(i, t1, tvel, pp, lastpp, t_force); //correct ...

						totalforce += sign*t_force;

					}
					else
					{
						//slipdir->slipvector2D !!!
						//totalforce += -sign*slipdir(i)*frictioncoeff*fabs(c_force)*(t1...t2); 
					}
				}

				GetMasterBody(j).GetdPosdqT(ploc,dpdq);

				for (int k=1; k <= f.Length(); k++)
				{
					f(k) -= -(dpdq(k,1)*totalforce.X()+dpdq(k,2)*totalforce.Y()+dpdq(k,3)*totalforce.Z());
				}
			}
		}

	}
};

Vector3D GeneralContact3D::GetRefPosD() const 
{
	Vector3D v;
	if (NSC() != 0)
	{
		v = GetSlaveNodePosD(1);
	}
	return Vector3D(v.X(),v.Y(),0);
}

void GeneralContact3D::DrawElement() 
{
	Vector3D p, pp;
	for (int i = 1; i <= NSC(); i++)
	{
		if (GetMBS()->GetIOption(112))
		{
			p = GetSlaveNodePosD(i);

			//draw contact node:
			mbs->SetColor(GetCol());
			//mbs->DrawSphere(GetSlaveBody(i).ToP3D(p),1*0.5*draw_dim.X(),4);

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

	for (int i = 1; i <= NMC(); i++)
	{
		//if (GetMBS()->GetIOption(112))
		if (masterisactive(i) == 0 || mathfunctions(masterisactive(i))->CheckCondition(GetMBS()->GetDrawTime()))
		{
			GeomElem(i)->DrawYourself();
		}
	}
};








