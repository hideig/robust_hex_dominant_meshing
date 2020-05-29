#include <TetraMeshTool/VertexInsertion.h>
#include <TetraMeshTool/StarBase.h>
//#include <TetraMeshTool/VertexNonsmooth.h>
#include <functional>
#include <iterator>
namespace VolumeMesh
{

//#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.26   //sin 15  0.2588
	#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.65  // sine 40  0.6428
//#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.45  // sine 25  0.4226

	/* fetch target tets which need vertex insertion */
	void VertexInsert::targetTetras(double &outminqual, int *trytetcnt, int *succcnt)
	{
		double minqual;/* quality of the worst tet in the mesh */
		TetraStack tetstack, outstack;      /* stack of tets  */
		double meshworstqual, origmeshworstqual;
		int numtets = 0;  /* the number of target tets*/

		/* figure out how many tets to take */
		//numtets = (int) (percentinsert * (double)tmesh->size_tetrahedron());

		/* calculate quality of each tet in mesh*/
		qualcalculator->meshquality(tmesh, qualitymetric, &minqual, &tetstack);

		/* sort the stack of tets from worst to best */
		sort(tetstack.begin(), tetstack.end(), comImproveTets);

		// get target tet counting
		for (unsigned int i = 0; i < tetstack.size(); i++)
		{
			if (tetstack[i].quality > VERTEX_INSERTION_QUALITY_THRESHOLD)
			{
				numtets = i;
				break;
			}
		}

		tetstack.erase(tetstack.begin()+numtets, tetstack.end());

		origmeshworstqual = meshworstqual = minqual;

		/* perform insertion pass on stack of tets */
		initImproveCommand();
		insertPass(tetstack, qualitymetric, outminqual, 0.2, succcnt);

		/* free the stack of tets */
		tetstack.clear();

		outminqual = qualcalculator->minmeshquality(tmesh, qualitymetric);

		if (trytetcnt != NULL)
			*trytetcnt = numtets;
	}


	/* attempt insertion on a stack of tets */
	bool VertexInsert::insertPass(TetraStack &tetstack, int qualmeasure, double &minqualafter, double okayqual, int *succcnt)
	{
		ImproveTetra * stacktet;
		int attempts = 0;
		int inserts = 0;
		bool success = false;
		TetraMesh::HedronVertexIter hv_iter;
		Point tetpoint[4];

		/* now go through the stack collecting information */
		while (!tetstack.empty())
		{
			/* pull the top tet off the stack */
			stacktet = &tetstack[0];
			tetstack.pop_front();

			/* make sure this tet still exists */
			if (!tmesh->is_valid(stacktet->handle)) 
				continue;

			/* because this tet might have been fixed up in a previous insertion,
			 * recompute its quality */
			hv_iter = tmesh->hedron_vertex_iter(stacktet->handle);
			tetpoint[0] = tmesh->point(hv_iter.handle());++hv_iter;
			tetpoint[1] = tmesh->point(hv_iter.handle());++hv_iter;
			tetpoint[2] = tmesh->point(hv_iter.handle());++hv_iter;
			tetpoint[3] = tmesh->point(hv_iter.handle());
			stacktet->quality = qualcalculator->tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualitymetric);

			/* if this tet is already above the demanded quality, don't attempt insertion */
			if (stacktet->quality >= okayqual && okayqual > 0.0)
				continue;

			/* attempt to insert a vertex */
			success = insertVertex(stacktet->handle);

			attempts++;
			inserts += success ? 1 : 0;

			success = false;
		}

		minqualafter = qualcalculator->minmeshquality(tmesh, qualitymetric);

		if (succcnt != NULL)
			*succcnt = inserts;

		if (inserts > 0) return true;

		return false;
	}

	bool VertexInsert::insertVertex(TetraMesh::HedronHandle tethandle, std::vector<TetraMesh::HedronHandle> *outtet)
	{
		PointHandle pnew;                      /* newly inserted point handle */
		std::vector<TetraMesh::HalfFaceHandle> newfaces;  /* faces of tets created after insert for seeding cavity drilling */
		std::vector<TetraMesh::HedronHandle> newTets;
		double cavityqual;       /* quality of the worst tet in the cavity */
		double worstdeletedqual;     /* quality of the worst tet deleted to form the cavity */
		Point p[4];
		Point barycenter(0,0,0);            /* barycenter of input tet */

		/* make sure this tet still exists */
		if (!tmesh->is_valid(tethandle)) return false;

		/* replace the bodyinsert */
		/* compute barycenter of this face */
		tetPoints(tmesh, tethandle, p);
		barycenter = p[1] + p[2] + p[3] + p[0];
		barycenter *= 0.25;

		/* allocate a new vertex */
		pnew = tmesh->add_point(barycenter, false);

		/* assign outputs */
		newfaces.clear();
		newTets.clear();
		newTets.push_back(tethandle);

		TetraMesh::HedronFaceIter hf_it;
		for (hf_it = tmesh->hedron_face_iter(tethandle); hf_it; ++hf_it)
		{
			newfaces.push_back(hf_it.handle());
		}
		//success = bodyinsert(tethandle, pnew, newfaces, newTets);
		/* end of replace */

		/* dig out optimal cavity */
		std::vector<CavityTet> outcavity;

		worstdeletedqual = HUGEFLOAT;
		cavityqual = HUGEFLOAT;

		/* build the cavity dag */
		buildcavitydag(pnew, newTets, newfaces, outcavity, true);

		if (maxcavitytet < outcavity.size())
			maxcavitytet = outcavity.size();

		/* build the cavity of maximum lexicographic quality */
		maxCavity(pnew, outcavity, newTets, &worstdeletedqual, &cavityqual);//#ifndef NO_TIMER

		/* did we succeed? */
		if (cavityqual > (worstdeletedqual + MININSERTIONIMPROVEMENT))
		{
			succvec.push_back(tethandle.idx());
			if (outtet != NULL)
			{
				outtet->clear();
				copy(newTets.begin(), newTets.end(), back_inserter(*outtet));
			}
			return true;
		}
		else
		{
			//rollback(&beforeid);
			tmesh->erase_point(pnew);
		}

		return false;
	}

	//bool VertexInsert::bodyinsert(TetraMesh::HedronHandle tetrahandle, PointHandle &newp, 
	//	                          std::vector<HalfFaceHandle> &newfaces, std::vector<TetraHandle> &newtets)
	//{
	//	PointHandle newtag;                  /* the tag of the new vertex */
	//	Point point[4];                     /* actual vertices of the tet */
	//	Point barycenter(0,0,0);            /* barycenter of input tet */
	//	std::vector<TetraHandle> tetVec;

	//	/* compute barycenter of this face */
	//	tetPoints(tmesh, tetrahandle, point);
	//	barycenter = point[1] + point[2] + point[3] + point[0];
	//	barycenter *= 0.25;

	//	/* allocate a new vertex */
	//	newtag = tmesh->add_point(barycenter, false);

	//	/* delete the old tet and add the three new ones */
	//	TMoperator->flip14(tetrahandle, barycenter, &tetVec);

	//	/* record this vertex insertion in the journal */
	//	pointvec.clear();thvec.clear();phvec.clear();vhvec.clear();
	//	TMoperator->getRecoverData(FLIP14, pointvec, thvec, phvec, vhvec);
	//	optjournal->pushJournal(FLIP14, pointvec, thvec, phvec, vhvec);

	//	/* check that we got it right */
	//	assert(tmesh->is_valid(tetVec[0]));
	//	assert(tmesh->is_valid(tetVec[1]));
	//	assert(tmesh->is_valid(tetVec[2]));
	//	assert(tmesh->is_valid(tetVec[3]));

	//	/* assign outputs */
	//	newp = newtag;
	//	newtets.push_back(tetVec[0]);
	//	newtets.push_back(tetVec[1]);
	//	newtets.push_back(tetVec[2]);
	//	newtets.push_back(tetVec[3]);

	//	TetraMesh::HedronFaceIter hf_iter;
	//	for (unsigned int i = 0; i < tetVec.size(); ++i)
	//	{
	//		hf_iter = tmesh->hedron_face_iter(tetVec[i]);
	//		for (; hf_iter; ++hf_iter)
	//			if (!tmesh->has_point(hf_iter.handle(), newtag))
	//			{
	//				newfaces.push_back(hf_iter.handle());
	//				break;
	//			}
	//	}


	//	/* its possible that insertion fails outright because the vertices may have moved
	//	a lot due to smoothing from their original positions. so, check if any of the new
	//	tets are inverted, and return false if that's the case. */
	//	for (int i=0; i<4; i++)
	//		if (tmesh->is_reverse(tetVec[i]))
	//			return false;

	//	return true;
	//}


	/************************************************************************/
	/*  Functions about Cavity                                              */
	/************************************************************************/

	/* find cavitytet by a tetrahedron handle*/
	unsigned int VertexInsert::findCavityTet(std::vector<CavityTet>::iterator begin, std::vector<CavityTet>::iterator end, TetraHandle handle)
	{
		unsigned int idx = 0;
		for (; begin != end; ++begin)
		{
			if (begin->handle == handle) break;
			++ idx;
		}
		return idx;
	}


	/* find cavityface int a cavityface container by a halfface handle */
	std::vector<CavityFace>::iterator VertexInsert::findCavityFace(std::vector<CavityFace>::iterator begin, std::vector<CavityFace>::iterator end, TetraMesh::HalfFaceHandle handle)
	{
		for (; begin != end; ++begin)
			if (begin->handle == handle)
				return begin;
		return end;
	}

	/* find tet adjacencies to a specified face 
	 * param[hf]: handle of the specified half face
	 * param[outin]: return the two adjacencies tets 
	                 outin[0]: tet contains the hf's opposite face
					 outin[1]: tet contains the face hf
	*/
	void VertexInsert::tetadjacencies(TetraMesh::HalfFaceHandle hf, TetraHandle *outin)
	{
		outin[1] = tmesh->handle_to_entity(hf).hedron_handle();
		if (tmesh->has_opposite_half_face(hf))
			outin[0] = tmesh->opposite_half_face(hf).hedron_handle();
		else
			outin[0] = GHOSTTET;
	}

	/* given a vertex, it's position, an initial cavitiy of tets and 
	 * a set of outward-oriented faces of that cavity, build a DAG
	 * representing the largest star-shaped cavity from the point of
	 * view of the inserted vertex */
	void VertexInsert::buildcavitydag(PointHandle pnew, std::vector<TetraHandle> initialC, std::vector<TetraMesh::HalfFaceHandle> initialFaces, std::vector<CavityTet> &outcavity, bool allowvertexdeletion)
	{
		std::vector<TetraMesh::HalfFaceHandle> F;    /* candidate face list */
		std::vector<TetraMesh::HalfFaceHandle> W;    /* wall face list */
		//std::vector<TetraHandle> C;                  /* cavity tet list */
		std::vector<TetraHandle> B;                  /* blocking tet list */
		TetraHandle currtet;                         /* current tet */
		TetraMesh::HalfFaceHandle currf;             /* current face */
		TetraHandle outin[2];
		TetraHandle outin2[2];
		Point tetpoint[4];
		int Wcount = 0;
		bool facing = true;
		TetraMesh::HalfFaceHandle otherfaces[3];
		int deepest = 0;
		std::vector<TetraMesh::HalfFaceHandle> nonwall;
		Point *tetpoints;
		TetraMesh::HedronFaceIter hf_iter;
		int depthlimit = cavdepthlimit;

		/* output cavity stuff */
		CavityTet cavtet;
		CavityFace cavface;
		CavityFace cavface2;
		unsigned int tetidx;

		/* fetch position of new vertex */
		const Point p = tmesh->point(pnew);

		/* initialize cavity tet list */
		for (unsigned int i = 0; i < initialC.size(); i++)
		{
			//C.push_back(initialC[i]);

			/* this tet is sure to be in the final cavity, add it */
			cavtet.handle = initialC[i];
			/* fake the qual because these tets were already split and made worse */
			cavtet.qual = 1.0;
			cavtet.depth = 0;
			cavtet.parents[0] = NOCAVITYTET;
			cavtet.parents[1] = NOCAVITYTET;
			cavtet.parents[2] = NOCAVITYTET;

			/* add cavity to the out put cavity vector */
			outcavity.push_back(cavtet);
		}

		/* initialize candidate face list */
		for (unsigned int i = 0; i < initialFaces.size(); i++)
		{
			F.push_back(initialFaces[i]);
		}

		/* now, as long as we have candidate faces */
		while (F.size() > 0)
		{
			assert(outcavity.size() <= MAXCAVITYTETS);
			assert(F.size() <= MAXCAVITYFACES);
			assert(W.size() <= MAXCAVITYFACES);
			//assert(C.size() <= MAXCAVITYTETS);
			assert(B.size() <= MAXCAVITYTETS);

			if (maxstackface < F.size())
				maxstackface = F.size();
			if (maxstackface < W.size())
				maxstackface = W.size();
			if (maxstacktet < B.size())
				maxstacktet = B.size();
			//if (maxstacktet < C.size())
			//	maxstacktet = C.size();

			/* pull a face out of F */
			currf = F[F.size()-1];
			F.pop_back();

			/* get t, the tet on the other side of this face 
			 * outin[0]: outward facing tet; outin[1]: inward facing tet*/
			tetadjacencies(currf, outin);

			/* the inward facing tet should already be in the output cavity.
			 * find it, and add this face as an outgoing face */
			tetidx = findCavityTet(outcavity.begin(), outcavity.end(), outin[1]);
			assert(tetidx != outcavity.size());

			/* compute the quality of the cavity tet with this face */
			Point fpoints[3];
			cavface.handle = currf;
			halfFacePoints(tmesh, currf, fpoints);
			cavface.qual = qualcalculator->tetquality(p, fpoints[0], fpoints[1], fpoints[2], qualitymetric);

			/* check to make sure it's not a ghost tet */
			if (outin[0] == GHOSTTET)
			{
				/* note that this face has no child, and assign it to its parent tet */
				cavface.child = NOCAVITYTET;
				if (outcavity[tetidx].outfaces.size() < 5)
					outcavity[tetidx].outfaces.push_back(cavface);

				/* add this face to the wall list */
				W.push_back(currf);
				continue;
			}

			// fetch the outward facing tet
			currtet = outin[0];

			// fetch the other faces of current tet
			int idx = 0;
			for (hf_iter = tmesh->hedron_face_iter(currtet); hf_iter; ++hf_iter)
			{
				/* except the opposite half face of the current face*/
				if (tmesh->opposite_half_face_handle(hf_iter.handle())==currf)
					continue;
				assert(idx < 3);
				otherfaces[idx ++] = hf_iter.handle();
			}

			/* fetch positions of vertices */
			tetPoints(tmesh, currtet, tetpoint);

			/* is t a cavity tet? */
			if (findCavityTet(outcavity.begin(), outcavity.end(), currtet) != outcavity.size())
			{
				/* we need to add this face to the parent tet, indicating 
				that it has no child because the tet on the other side
				doesn't depend in it's removal to exist */
				cavface.child = NOCAVITYTET;
				if (outcavity[tetidx].outfaces.size() < 5)
					outcavity[tetidx].outfaces.push_back(cavface);

				/* yes, do nothing */
				continue;
			}

			/* is t a blocking tet? */
			if (find(B.begin(), B.end(), currtet) != B.end())
			{
				/* if there is one other wall face of this tet, and the other two faces
				are visible from v, we can add this tet to the cavity. */
				Wcount = 1;
				facing = true;

				TetraMesh::HedronFaceIter hf_iter;
				for (hf_iter = tmesh->hedron_face_iter(currtet); hf_iter; ++hf_iter)
				{
					/* except the current face */
					if (hf_iter.handle() == tmesh->opposite_half_face_handle(currf))
						continue;

					/* is this face already marked as a wall ? 
					 * testing the opposite face which may be visited before */
					if (find(W.begin(), W.end(), tmesh->opposite_half_face_handle(hf_iter.handle())) != W.end())
					{
						Wcount++;
					}
					else
					/* it's not a wall... is it oriented toward v? */
					{
						Point fp[3];
						halfFacePoints(tmesh, hf_iter.handle(), fp);
						if (orient(p, fp[0], fp[1], fp[2]) <= MINFACING)
						{
							facing = false;
						}
					}
				}

				/* only allow tets with three parents if we are allowing vertex deletion */
				if ((Wcount == 2 || (Wcount == 3 && allowvertexdeletion)) && facing)
				{
					/* this tet can be added to the cavity */

					/* remove it from B */
					B.erase(find(B.begin(), B.end(), currtet));
					/* add it to C */
					//C.push_back(currtet);

					/* add this tet to the output cavity */
					cavtet.handle = currtet;

					/* compute it's original quality */
					tetpoints = new Point[4];
					tetPoints(tmesh, currtet, tetpoints);
					cavtet.qual = qualcalculator->tetquality( tetpoints[0], tetpoints[1], tetpoints[2], tetpoints[3], qualitymetric);
					delete [] tetpoints;

					/* we know one parent must be the one we found above */
					cavtet.parents[0] = outcavity[tetidx].handle;
					/* the depth is one more than the parent depth */
					cavtet.depth = outcavity[tetidx].depth + 1;
					/* if this is a new deepest, remember it */
					if (cavtet.depth > deepest) deepest = cavtet.depth;
					outcavity.push_back(cavtet);

					/* add this face to the parent tet with the correct child */
					//cavface.child = cavtet;
					if (outcavity[tetidx].outfaces.size() < 5)
						outcavity[tetidx].outfaces.push_back(cavface);

					/* remove any faces that were in W, add others to F. Handle output
					tet faces that need to be added */
					nonwall.clear();

					/* first, handle all wall face so we can set the correct depth */
					TetraMesh::HedronFaceIter hf_iter;
					for (hf_iter = tmesh->hedron_face_iter(currtet); hf_iter; ++ hf_iter)
					{
						/* except the current face */
						if (hf_iter.handle() == tmesh->opposite_half_face_handle(currf))
							continue;

						/* is this already a wall face? */
						std::vector<TetraMesh::HalfFaceHandle>::iterator w_it = find(W.begin(), W.end(), tmesh->opposite_half_face_handle(hf_iter.handle()));
						if (w_it != W.end())
						{
							W.erase(w_it);

							/* because this face was previously a wall face,
							it has some cavity tet that it belongs to. find
							this tet in the output cavity and set it's child face */
							unsigned int ctidx;
							ctidx = findCavityTet(outcavity.begin(), outcavity.end(), tmesh->handle_to_entity(tmesh->opposite_half_face_handle(hf_iter.handle())).hedron_handle());
							assert(ctidx != outcavity.size());

							/* add this face to the parent tet's outgoing faces */
							cavface.handle = tmesh->opposite_half_face_handle(hf_iter.handle());
							halfFacePoints(tmesh, cavface.handle, fpoints);
							cavface.qual = qualcalculator->tetquality( tmesh->point(pnew), fpoints[0], fpoints[2], fpoints[1], qualitymetric);
							cavface.child = currtet;

							/* make sure that this face is already in this tet */
							if (outcavity[ctidx].outfaces.size() < 5)
								outcavity[ctidx].outfaces.push_back(cavface);

							/* assign the parent tet as the second parent of the new cavity tet */
							if (outcavity[outcavity.size() -1].parents[1] == NOCAVITYTET)
							{
								outcavity[outcavity.size() -1].parents[1] = outcavity[ctidx].handle;
							}
							else
							{
								assert(outcavity[outcavity.size() -1].parents[2] == NOCAVITYTET);
								outcavity[outcavity.size() -1].parents[2] = outcavity[ctidx].handle;
							}

							/* if this parent has a lesser depth value, update new tet's depth to be the lesser */
							if (outcavity[ctidx].depth < outcavity[outcavity.size() -1].depth)
							{
								outcavity[outcavity.size() -1].depth = outcavity[ctidx].depth;
							}
						}
						else
						{
							/* record this non-wall face for potential addition to F later */
							nonwall.push_back(hf_iter.handle());
						}
					}

					for (unsigned int i = 0; i < nonwall.size(); ++i)
					{
						/* this is a newly-uncovered face. there could be more tets behind it, so
						we should add it to F, if the current tet's depth isn't more than the max */
						if (outcavity[outcavity.size() -1].depth < depthlimit)
						{
							F.push_back(nonwall[i]);
						}
						/* we should artificially make this a wall face so the cavity doesn't get deeper */
						else
						{
							/* construct output face */
							Point fp[3];
							cavface2.handle = nonwall[i];
							halfFacePoints(tmesh, nonwall[i], fp);
							cavface2.qual = qualcalculator->tetquality(p, fp[0], fp[2], fp[1], qualitymetric);
							cavface2.child = NOCAVITYTET;
							//assert(cavface2.qual > 0);

							/* add it to parent tet */
							if ((outcavity[outcavity.size()-1].outfaces).size() < 5)
								(outcavity[outcavity.size()-1].outfaces).push_back(cavface2);

							W.push_back(nonwall[i]);
						}
					}
				}
				else
				{
					/* note that this face has no child, and assign it to its parent tet */
					cavface.child = NOCAVITYTET;
					if (outcavity[tetidx].outfaces.size() < 5)
						outcavity[tetidx].outfaces.push_back(cavface);

					/* add f to W, it borders a blocking tet */
					W.push_back(currf);
				}
				continue;
			}

			/* t is neither a blocking tet nor a cavity tet */
			/* check to see if the three other faces of the tet are facing v */
			bool check = true;
			for (int i = 0; i < 3; ++i)
			{
				/* fetch the face points */
				Point fp[3];
				halfFacePoints(tmesh, otherfaces[i], fp);
				/* the order of points should be reversed */
				if (orient(p, fp[2], fp[1], fp[0]) != 1)
				{
					check = false;
					break;
				}
			}
			if (check)
			{
				/* yes! we can add this tet to the cavity */
				//C.push_back(currtet);

				/* add this tet to the output cavity */
				cavtet.handle = currtet;
				/* compute it's original quality */
				tetpoints = new Point[4];
				tetPoints(tmesh, currtet, tetpoints);
				cavtet.qual = qualcalculator->tetquality(tetpoints[0], tetpoints[1], tetpoints[2], tetpoints[3], qualitymetric);
				delete [] tetpoints;
				/* it's parent must be the parent above */
				cavtet.parents[0] = outcavity[tetidx].handle;
				/* depth is one deeper than parent */
				cavtet.depth = outcavity[tetidx].depth + 1;
				/* if this is a new deepest, note it */
				if (cavtet.depth > deepest) deepest = cavtet.depth;
				outcavity.push_back(cavtet);

				/* note the current face's child in the parent tet */
				cavface.child = currtet;
				if (outcavity[tetidx].outfaces.size() < 5)
					outcavity[tetidx].outfaces.push_back(cavface);

				/* add t's three (outward oriented) faces to F, if the current tet isn't too deep */
				if (cavtet.depth < depthlimit)
				{
					F.push_back(otherfaces[0]);
					F.push_back(otherfaces[1]);
					F.push_back(otherfaces[2]);
				}
				else
				{
					/* construct output face */
					cavface2.child = NOCAVITYTET;

					Point fp[3];
					cavface2.handle = otherfaces[0];
					halfFacePoints(tmesh, otherfaces[0], fp);
					cavface2.qual = qualcalculator->tetquality(p, fp[0], fp[1], fp[2], qualitymetric);
					/* add it to parent tet */
					(outcavity[outcavity.size()-1].outfaces).push_back(cavface2);

					cavface2.handle = otherfaces[1];
					halfFacePoints(tmesh, otherfaces[1], fp);
					cavface2.qual = qualcalculator->tetquality(p, fp[0], fp[1], fp[2], qualitymetric);
					/* add it to parent tet */
					(outcavity[outcavity.size()-1].outfaces).push_back(cavface2);

					cavface2.handle = otherfaces[2];
					halfFacePoints(tmesh, otherfaces[2], fp);
					cavface2.qual = qualcalculator->tetquality(p, fp[0], fp[1], fp[2], qualitymetric);
					/* add it to parent tet */
					(outcavity[outcavity.size()-1].outfaces).push_back(cavface2);

					W.push_back(otherfaces[0]);
					W.push_back(otherfaces[1]);
					W.push_back(otherfaces[2]);
				}
			}
			else
			{
				/* this is a blocking tet, add it to B */
				B.push_back(currtet);

				/* note the current face in the parent tet */
				cavface.child = NOCAVITYTET;
				if ((outcavity[tetidx].outfaces).size() < 5)
				{
					(outcavity[tetidx].outfaces).push_back(cavface);
				}

				/* add the current face to the wall face list */
				W.push_back(currf);
			}
		}

		/* record the maximum depth */
		cavdeep = deepest;
	}


	/* add an outgoing face to a cavity tet */
	void VertexInsert::addCavityTetFace(CavityTet &tet, CavityFace &cface)
	{
		int parentcount = 0;

		/* check if this face already exists (i.e., we're updating it) */
		std::vector<CavityFace>::iterator outf_iter;
		outf_iter = findCavityFace(tet.outfaces.begin(), tet.outfaces.end(), cface.handle);
		
		/* we can't add a new face if there are already 3 */
		if (outf_iter == tet.outfaces.end())
		{
			/* make sure that the sum of the parents and outbound faces is 4 */
			parentcount = numParents(tet);
			/* allow more than 3 outfaces for parentless tets */
			if (parentcount != 0)
			{
				assert(parentcount + tet.outfaces.size() < 4);
			}
			assert(tet.outfaces.size() <= MAXOUTFACES);
			/* copy in the face */
			tet.outfaces.push_back(cface);
		}
	}

	/* return the number of parents a tet has */
	int VertexInsert::numParents(CavityTet &tet)
	{
		int numparents = 0;
		if (tet.parents[0] != NOCAVITYTET) numparents++;
		if (tet.parents[1] != NOCAVITYTET) numparents++;
		if (tet.parents[2] != NOCAVITYTET) numparents++;
		return numparents;
	}

	/* ALGORITHM 3 + CONSIDER DELETED TETS*/
	/* given the DAG representing the largest possible cavity,
	find the subcavity that has the lexicographically maximum
	quality, then dig it out and connect the exposed faces to 
	the new vertex, returning the new tets as well as the worst
	quality in this cavity */
	void VertexInsert::maxCavity(PointHandle pnew, std::vector<CavityTet> &cavity, std::vector<TetraHandle> &outtets, 
		                         double *worstdeleted, double *worstincavity)
	{
		bool foundparent = false;                /* did we find face parent? */
		std::vector<CavityEdgeorTet> edges;      /* list of all edges in the DAG */
		CavityTet t;                   /* the virtual node t, at the end of array */
		int parentlabel, childlabel;             /* the groups that contain parent and child of an edge */
		double depthfactor;
		double qual;
		int deepest = 0;
		unsigned int parenttetidx, childtetidx;
		CavityEdgeorTet tmpedge;
		std::vector<TetraMesh::HalfFaceHandle> outputfaces;
		std::vector<TetraHandle> erasetetras;  /* tetras which are going to be erased */
		unsigned int i;

		/* initialize t as a ghost tet */
		t.handle = TetraHandle(-1);

		/* now, proceed through the DAG, recording all edges, and deleting child
		information to remove the edges, producing D' */
		for (i = 0; i < cavity.size(); i++)
		{
			/* initialize this tet's label. set s = cavity, others no label */
			if (i == 0)
			{
				cavity[i].label = CAVLABEL;
			}
			else
			{
				cavity[i].label = NOLABEL;
			}

			/* compute how much to exagerrate quality because of depth */
			assert(cavity[i].depth >= 0);
			if (cavity[i].depth < DEPTHTABLESIZE)
			{
				depthfactor = cavitydepthtable[cavity[i].depth];
			}
			else
			{
				depthfactor = cavitydepthtable[DEPTHTABLESIZE - 1];
			}
			assert(depthfactor >= 1.0);

			/* for each outgoing face */
			for (unsigned int j = 0; j < cavity[i].outfaces.size(); j++)
			{
				tmpedge.label = EDGELABEL;
				tmpedge.parent = cavity[i].handle;
				tmpedge.qual = cavity[i].outfaces[j].qual * depthfactor;
				/* save which outgoing face this child was */
				tmpedge.childnum = j;

				if (cavity[i].outfaces[j].child == NOCAVITYTET)
				{
					/* if it has no child tet, make up a virtual edge to t */
					tmpedge.child = GHOSTTET;
				}
				else
				{
					/* set edge child to be actual child */
					tmpedge.child = cavity[i].outfaces[j].child;
				}
				edges.push_back(tmpedge);
				assert(edges.size() < MAXCAVITYDAGEDGES);

				/* initialize H by setting no edges to be in it */
				cavity[i].outfaces[j].inH = false;
			}
		}

		if (maxcavityedge < edges.size())
			maxcavityedge = edges.size();
		/* now, sort the edges in order of ascending quality */
		sort(edges.begin(), edges.end(), comCavityEdge_less);

		/* go through each edge, adding it to D' if it doesn't
		connect s to t */

		std::vector<int> labeltet; // for testing
		std::vector<int> antilabeltet;
		std::vector<int> labeltetloop;
		std::vector<int> antilabeltetloop;
		std::vector<CavityTet> cavityloop;
		std::vector<TetraHandle> erasetetrasloop;
		copy(cavity.begin(), cavity.end(), back_inserter(cavityloop));
		for (i = 0; i < edges.size(); i++)
		{
			/* find parent cavity */
			parenttetidx = findCavityTet(cavity.begin(), cavity.end(), edges[i].parent);

			/* check parent's label */
			parentlabel = cavity[parenttetidx].label;
			/* check child's label */
			childtetidx = findCavityTet(cavity.begin(), cavity.end(), edges[i].child);
			if (edges[i].child == GHOSTTET)
				/* this child is labeled automatically as anti-cavity */
				childlabel = ANTICAVLABEL;
			else
				childlabel = cavity[childtetidx].label;

			/* if the parent is in the cavity */
			if (parentlabel == CAVLABEL)
			{
				/* and the child is in the anti-cavity */
				if (childlabel == ANTICAVLABEL)
				{
					/* record output face from parent to child, record the half face of child */
					outputfaces.push_back(cavity[parenttetidx].outfaces[edges[i].childnum].handle);
				}
				/* otherwise, if the child isn't labeled */
				else
				{
					if (childlabel == NOLABEL)
					{
						//labeltet.clear();
						//cavityLabel(cavity, childtetidx, &labeltet);
						//labeltetloop.clear();
						//cavityLabel_loop(cavityloop, childtetidx, labeltetloop);

						labeltet.clear();
						cavityLabel(cavityloop, childtetidx, &labeltet);
						labeltetloop.clear();
						cavityLabel_loop(cavity, childtetidx, labeltetloop);
					}
				}
			}
			/* parent isn't labeled cavity. */
			else
			{
				/* is the parent wholly unlabeled ? */
				if (parentlabel == NOLABEL)
				{
					/* is the child labeled anti-cavity ? */
					if (childlabel == ANTICAVLABEL)
					{
						//antilabeltet.clear();
						//antiCavityLabel(cavity, parenttetidx, &antilabeltet);
						//antilabeltetloop.clear();
						//antiCavityLabel_loop(cavityloop, parenttetidx, antilabeltetloop);

						antilabeltet.clear();
						antiCavityLabel(cavityloop, parenttetidx, &antilabeltet);
						antilabeltetloop.clear();
						antiCavityLabel_loop(cavity, parenttetidx, antilabeltetloop);
					}
					else
					/* neither the parent nor the child is labeled */
					{
						/* add the edge from parent to child to H */
						cavity[parenttetidx].outfaces[edges[i].childnum].inH = true;
					}
				}
			}
		}

		/* keep track of what the deepest tet in the final cavity was */
		deepest = 0;

		/* delete all tets labeled as cavity */
		Point p[4];
		*worstdeleted = HUGEFLOAT;
		for (i = 0; i < cavity.size(); i++)
		{
			if (cavity[i].label == CAVLABEL)
			{
				tetPoints(tmesh, cavity[i].handle, p);
				cavity[i].qual = qualcalculator->tetquality(p[0], p[1], p[2], p[3], qualitymetric);

				/* is this the worst quality tet we're deleting? */
				if (cavity[i].qual < *worstdeleted)
				{
					*worstdeleted = cavity[i].qual;
				}

				/* is this the deepest tet we've encountered? */
				if (cavity[i].depth > deepest)
				{
					deepest = cavity[i].depth;
				}
				erasetetras.push_back(cavity[i].handle);
			}

			if (cavityloop[i].label == CAVLABEL)
			{
				tetPoints(tmesh, cavityloop[i].handle, p);
				cavityloop[i].qual = qualcalculator->tetquality(p[0], p[1], p[2], p[3], qualitymetric);

				/* is this the worst quality tet we're deleting? */
				if (cavityloop[i].qual < *worstdeleted)
				{
					*worstdeleted = cavityloop[i].qual;
				}

				/* is this the deepest tet we've encountered? */
				if (cavityloop[i].depth > deepest)
				{
					deepest = cavityloop[i].depth;
				}
				erasetetrasloop.push_back(cavityloop[i].handle);
			}
		}

		/* record depth of deepest tet */
		cavdeep = deepest;

		if (maxcavityface < outputfaces.size())
			maxcavityface = outputfaces.size();

		double newworstquality = HUGEFLOAT;
		double tempquality;
		/* now build the output cavity using the output faces */
		std::vector<Point> facespoints;
		outtets.clear();
		for (unsigned int i = 0; i < outputfaces.size(); i++)
		{
			assert(tmesh->is_valid(outputfaces[i]));
			Point fp[3];
			halfFacePoints(tmesh, outputfaces[i], fp);
			facespoints.push_back(fp[0]);
			facespoints.push_back(fp[1]);
			facespoints.push_back(fp[2]);	
		}


		// smooth insert point
		Point newpoint;
		newpoint = tmesh->point(pnew);
		smoothInsertVertex(newpoint, facespoints, qualitymetric);
		tmesh->set_point(pnew, newpoint);

		// calculate new local quality
		for (unsigned int i = 0; i < outputfaces.size(); i++)
		{
			tempquality = qualcalculator->tetquality(newpoint, facespoints[i*3], facespoints[i*3+1], facespoints[i*3+2], qualitymetric);
			if (newworstquality > tempquality)
				newworstquality = tempquality;
		}
	
		*worstincavity = newworstquality;
		if (newworstquality - *worstdeleted < 1.0*10e-6)
		{		
			return ;
		}

		/* erase original tetras */
		for (unsigned int i = 0; i < erasetetras.size();++i)
		{
			tmesh->erase_tetrahedron(erasetetras[i]);
		}
		/* end of tet erasing */

		/* add new tetras */
		for (unsigned int i = 0; i < facespoints.size(); i += 3)
		{
			outtets.push_back(tmesh->add_tetrahedron(tmesh->point(pnew), facespoints[i], facespoints[i+1], facespoints[i+2]));
		}
		/* end of tet adding */
	}

	// smooth the inserting vertex with the cavity faces information
	bool VertexInsert::smoothInsertVertex(Point &p, std::vector<Point> &faces, int qualmetric)
	{
		Point newp(0,0,0);
		double qualbefore, qualafter;
		double tempqual;
		//double point[3];
		//double **tetp;
		int tetcnt;
		//VolumeMesh::Nonsmoother nonsmoother;
		std::set<Point> neipoint;
		std::set<Point>::iterator neipoint_it;
		qualbefore = qualafter = 1.0;

		// calculate local quality before smoothing
		for (unsigned int i = 0; i < faces.size(); i = i+3)
		{
			neipoint.insert(faces[i]);
			neipoint.insert(faces[i+1]);
			neipoint.insert(faces[i+2]);
			tempqual = qualcalculator->tetquality(p, faces[i], faces[i+1], faces[i+2], qualmetric);
			if (qualbefore > tempqual)
				qualbefore = tempqual;
		}

		/* smooth the point */
		for (neipoint_it = neipoint.begin(); neipoint_it != neipoint.end(); ++neipoint_it)
		{
			newp += *neipoint_it;
		}
		newp /= neipoint.size();
		// fetch smoothing information
		//tetcnt = faces.size()/3;
		//tetp = new double *[tetcnt];
		//for (int i = 0; i < tetcnt; i++)
		//{
		//	tetp[i] = new double[9];
		//	for (int j = 0; j < 3; j++)
		//	{
		//		tetp[i][j*3]   = faces[i*3+j][0];
		//		tetp[i][j*3+1] = faces[i*3+j][1];
		//		tetp[i][j*3+2] = faces[i*3+j][2];
		//	}
		//}
		//point[0] = p[0]; point[1] = p[1]; point[2] = p[2];

		// do vertice smoothing
		//nonsmoother.setTmesh(tmesh);
		//nonsmoother.setQualityKind(MINSINE2);
		//nonsmoother.NonsmoothVertice(point, tetp, tetcnt, point);
		/* end of smoothing */

		// calculate local quality after smoothing
		//newp[0] = point[0];newp[1] = point[1];newp[2] = point[2];
		for (unsigned int i = 0; i < faces.size(); i = i+3)
		{
			tempqual = qualcalculator->tetquality(newp, faces[i], faces[i+1], faces[i+2], qualmetric);
			if (qualafter > tempqual)
				qualafter = tempqual;
		}

		if (qualafter > qualbefore)
		{
			p = newp;
			return true;
		}
		return false;
	}

	/* recursively label parents and children as cavity tets */
	void VertexInsert::cavityLabel(std::vector<CavityTet> &cavity, unsigned int ctetidx, std::vector<int> *labeltets)
	{
		unsigned int parenttetidx;
		unsigned int childtetidx;
		/* this tet shouldn't yet be labeled */
		if (cavity[ctetidx].label != NOLABEL)
			return;
		assert(cavity[ctetidx].label == NOLABEL);

		/* label this tet as in the cavity */
		cavity[ctetidx].label = CAVLABEL;

		if (labeltets != NULL)
			labeltets->push_back(ctetidx);

		/* go through all parents in the original graph */
		for (int i = 0; i < 3; ++i)
		{
			if (cavity[ctetidx].parents[i] != NOCAVITYTET)
			{
				/* if this parent is unlabeled, label it */
				parenttetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].parents[i]);
				if (parenttetidx != cavity.size() && cavity[parenttetidx].label == NOLABEL)
				{
					cavityLabel(cavity, parenttetidx);
				}
			}
		}

		/* go through all children in H */
		for (unsigned int i = 0;  i < cavity[ctetidx].outfaces.size(); ++i)
		{
			/* check if this edge is in H */
			if (cavity[ctetidx].outfaces[i].inH == true)
			{
				/* this can't be an edge leading to t... we should never add those */
				assert(cavity[ctetidx].outfaces[i].child != NOCAVITYTET);

				childtetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].outfaces[i].child);
				if (childtetidx != cavity.size() && cavity[childtetidx].label == NOLABEL)
				{
					cavityLabel(cavity, childtetidx);
				}
			}
		}
	}


	/* recursively label parents and children as cavity tets */
	void VertexInsert::cavityLabel_loop(std::vector<CavityTet> &cavity, unsigned int ctetidx_, std::vector<int> &labeltets)
	{
		unsigned int parenttetidx;
		unsigned int childtetidx;
		unsigned int tetstack[MAXCAVITYTETS];
		unsigned int tetstacksize = 0;
		unsigned int ctetidx;
		int i;

		/* this tet shouldn't yet be labeled */
		if (cavity[ctetidx_].label != NOLABEL)
			return;
		assert(cavity[ctetidx_].label == NOLABEL);

		tetstack[tetstacksize++] = ctetidx_;

		while (tetstacksize > 0)
		{
			if (maxstacktet < tetstacksize)
				maxstacktet = tetstacksize;

			ctetidx = tetstack[tetstacksize-1];
			-- tetstacksize;

			/* label this tet as in the cavity */
			labeltets.push_back(ctetidx);
			if (cavity[ctetidx_].label != NOLABEL)
				continue;
			cavity[ctetidx].label = CAVLABEL;

			/* go through all parents in the original graph */
			for (i = 0; i < 3; ++i)
			{
				if (cavity[ctetidx].parents[i] != NOCAVITYTET)
				{
					/* if this parent is unlabeled, label it */
					parenttetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].parents[i]);
					if (parenttetidx != cavity.size() && cavity[parenttetidx].label == NOLABEL)
					{
						tetstack[tetstacksize++] = parenttetidx;
					}
				}
			}


			/* go through all children in H */
			for (i = 0;  i < cavity[ctetidx].outfaces.size(); ++i)
			{
				/* check if this edge is in H */
				if (cavity[ctetidx].outfaces[i].inH == true)
				{
					/* this can't be an edge leading to t... we should never add those */
					assert(cavity[ctetidx].outfaces[i].child != NOCAVITYTET);

					childtetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].outfaces[i].child);
					if (childtetidx != cavity.size() && cavity[childtetidx].label == NOLABEL)
					{
						tetstack[tetstacksize++] = childtetidx;
					}
				}
			}
		}
	}


	/* recursively label parents and children as anti-cavity tets */
	void VertexInsert::antiCavityLabel(std::vector<CavityTet> &cavity, unsigned int ctetidx, std::vector<int> *labeltets)
	{
		unsigned int parenttetidx;
		unsigned int childtetidx;
		int parent,edgetochild;

		/* this tet shouldn't yet be labeled */
		if (cavity[ctetidx].label != NOLABEL)
			return;
		assert(cavity[ctetidx].label == NOLABEL);

		/* label this tet as in the anticavity */
		cavity[ctetidx].label = ANTICAVLABEL;
		if (labeltets != NULL)
			labeltets->push_back(ctetidx);

		/* go through all parents in H */
		for (int i = 0; i < 3; ++i)
		{
			parent = cavity[ctetidx].parents[i];

			if (parent != NOCAVITYTET)
			{
				/* is this parent unlabeled ? */
				parenttetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].parents[i]);
				if (parenttetidx == cavity.size() && cavity[parenttetidx].label != NOLABEL)
				{
					continue;
				}

				/* find the edge from this parent down to the child */
				edgetochild = -1;
				for (unsigned int j = 0; j < cavity[parenttetidx].outfaces.size(); j++)
				{
					if (cavity[parenttetidx].outfaces[j].child == cavity[ctetidx].handle)
					{
						edgetochild = j;
					}
				}
				/* make sure we found the edge */
				assert(edgetochild != -1);

				/* is this edge in H? */
				if (cavity[parenttetidx].outfaces[edgetochild].inH == true)
				{
					antiCavityLabel(cavity, parenttetidx);
				}
			}
		}

		/* go through all children in original graph G */
		for (unsigned i = 0; i < cavity[ctetidx].outfaces.size(); i++)
		{
			/* if the child is t, it's the end and is already labeled. move on */
			if (cavity[ctetidx].outfaces[i].child == NOCAVITYTET)
			{
				continue;
			}

			childtetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].outfaces[i].child);
			if (childtetidx != cavity.size() && cavity[childtetidx].label == NOLABEL)
			{
				antiCavityLabel(cavity, childtetidx);
			}
		}
	}


	/* recursively label parents and children as anti-cavity tets */
	void VertexInsert::antiCavityLabel_loop(std::vector<CavityTet> &cavity, unsigned int ctetidx_, std::vector<int> &labeltets)
	{
		unsigned int parenttetidx;
		unsigned int childtetidx;
		int parent,edgetochild;
		unsigned int tetstack[MAXCAVITYTETS];
		unsigned int tetstacksize = 0;
		unsigned int ctetidx;
		int i;

		/* this tet shouldn't yet be labeled */
		if (cavity[ctetidx_].label != NOLABEL)
			return;
		assert(cavity[ctetidx_].label == NOLABEL);

		tetstack[tetstacksize++] = ctetidx_;

		while (tetstacksize > 0)
		{
			if (maxstacktet < tetstacksize)
				maxstacktet = tetstacksize;

			ctetidx = tetstack[tetstacksize-1];
			--tetstacksize;

			/* label this tet as in the anticavity */
			if (cavity[ctetidx_].label != NOLABEL)
				continue;
			cavity[ctetidx].label = ANTICAVLABEL;
			labeltets.push_back(ctetidx);

			/* go through all parents in H */
			for (i = 0; i < 3; ++i)
			{
				parent = cavity[ctetidx].parents[i];

				if (parent != NOCAVITYTET)
				{
					/* is this parent unlabeled ? */
					parenttetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].parents[i]);
					if (parenttetidx == cavity.size() && cavity[parenttetidx].label != NOLABEL)
					{
						continue;
					}

					/* find the edge from this parent down to the child */
					edgetochild = -1;
					for (unsigned int j = 0; j < cavity[parenttetidx].outfaces.size(); j++)
					{
						if (cavity[parenttetidx].outfaces[j].child == cavity[ctetidx].handle)
						{
							edgetochild = j;
						}
					}
					/* make sure we found the edge */
					assert(edgetochild != -1);

					/* is this edge in H? */
					if (cavity[parenttetidx].outfaces[edgetochild].inH == true)
					{
						tetstack[tetstacksize++] = parenttetidx;
					}
				}
			}

			/* go through all children in original graph G */
			for (unsigned i = 0; i < cavity[ctetidx].outfaces.size(); i++)
			{
				/* if the child is t, it's the end and is already labeled. move on */
				if (cavity[ctetidx].outfaces[i].child == NOCAVITYTET)
				{
					continue;
				}

				childtetidx = findCavityTet(cavity.begin(), cavity.end(), cavity[ctetidx].outfaces[i].child);
				if (childtetidx != cavity.size() && cavity[childtetidx].label == NOLABEL)
				{
					tetstack[tetstacksize++] = childtetidx;
				}
			}
		}
	}

	/************************************************************************/
	/* Convenience functions                                                */
	/************************************************************************/

	/* A predicate function for cavityedge comparision */
	bool VertexInsert::comCavityEdge_less(const CavityEdgeorTet &e1, const CavityEdgeorTet &e2)
	{
		return e1.qual < e2.qual;
	}
	/* A predicate function for cavitytet comparison*/
	bool VertexInsert::comCavityTet(const CavityTet &a, const CavityTet &b)
	{
		return a.handle == b.handle;
	}

	bool VertexInsert::comCavityTetHandle(const CavityTet &a, const TetraHandle & handle)
	{
		return a.handle == handle;
	}
	/* A predicate function for ImproveTetra comparison*/
	bool VertexInsert::comImproveTets(const ImproveTetra &a, const ImproveTetra &b)
	{
		return a.quality < b.quality;
	}

	/* get four points of a tetrahedron indexed by a TetraHandle*/
	bool VertexInsert::tetPoints(TetraMesh *mesh, TetraHandle handle, Point *p)
	{
		if (p == NULL)
			p = new Point[4];
		if (handle.idx()<0 || handle.idx() >= (int)mesh->size_tetrahedron() )
			return false;
		TetraMesh::HedronVertexIter hv_iter;
		hv_iter = mesh->hedron_vertex_iter(handle);
		p[0] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[1] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[2] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[3] = mesh->point(hv_iter.handle());
		return true;
	}

	/* get three points of a half face indexed by a TetraHandle*/
	bool VertexInsert::halfFacePoints(TetraMesh *tmesh, TetraMesh::HalfFaceHandle hf, Point *p)
	{
		if (p == NULL)
			p = new Point [3];
		if (hf.idx()<0 || hf.idx() >= (int)tmesh->size_halfface())
			return false;
		TetraMesh::HalfFaceVertexIter fv_iter;
		fv_iter = tmesh->half_face_vertex_iter(hf);
		p[0] = tmesh->point(fv_iter.handle()); ++fv_iter;
		p[1] = tmesh->point(fv_iter.handle()); ++fv_iter;
		p[2] = tmesh->point(fv_iter.handle());
		return true;
	}

	/* get three points of a half face indexed by a TetraHandle*/
	bool VertexInsert::halfFacePointHandle(TetraMesh *tmesh, TetraMesh::HalfFaceHandle hf, PointHandle *p)
	{
		if (p == NULL)
			p = new PointHandle [3];
		if (hf.idx()<0 || hf.idx() >= (int)tmesh->size_halfface())
			return false;
		TetraMesh::HalfFaceVertexIter fv_iter;
		fv_iter = tmesh->half_face_vertex_iter(hf);
		p[0] = tmesh->point_handle(fv_iter.handle()); ++fv_iter;
		p[1] = tmesh->point_handle(fv_iter.handle()); ++fv_iter;
		p[2] = tmesh->point_handle(fv_iter.handle());
		return true;
	}

	void VertexInsert::initImproveCommand()
	{
		/* insertion options */
		insertthreshold = 0.035;     /* percent worst tets */
		cavityconsiderdeleted = 0;   /* consider enlarging cavity for deleted tets? */
		cavdepthlimit = 6;           /* only allow initial cavity to includes tets this deep */
	}
}