#include <TetraMeshTool/VertexSmoothing.h>
#include <TetraMeshTool/top.h>
#include <stack>
#include <TetraMeshTool/StarBase.h>
#include <iterator>
namespace VolumeMesh
{
	void VertexSmoother::initialize(int do_nonsmooth_, int fixedsmooth_)
	{
		/* set smoothing data */
		sinewarpfactor = 0.75;
		quadricoffset = 0.8;
		quadricscale = 300.0;

		/*smoothing options*/
		do_nonsmooth = do_nonsmooth_;
		fixedsmooth = fixedsmooth_;
		usequadrics = 1;

		segment();
	}

	/* Do the segmentation of the model
	 * Mark the half faces with the segment they are in
	*/
	void VertexSmoother::segment()
	{
		std::stack<TetraMesh::HalfFaceHandle> grow;
		std::stack<TetraMesh::HalfFaceHandle> seedstack;
		TetraMesh::HalfFaceHandle seed,curr,newface;
		TetraMesh::HedronIter h_it;
		TetraMesh::HedronFaceIter hf_it;
		TetraMesh::FaceHalfedgeIter fe_it;
		TetraMesh::HalfEdgeHandle mate;
		bool newregion = true;
		regionnum = 1;

		// fetch the first boundary half face
		for (h_it = tmesh->hedrons_begin(); h_it != tmesh->hedrons_end(); ++h_it)
		{
			if (!tmesh->is_valid(h_it.handle()))
				continue;
			for (hf_it = tmesh->hedron_face_iter(h_it.handle()); hf_it; ++hf_it)
			{
				if (tmesh->is_boundary(hf_it.handle()))
				{
					seed = hf_it.handle();
					break;
				}
			}
			if (seed != TetraMesh::HalfFaceHandle(-1))
				break;
		}
		seedstack.push(seed);

		// start a segmentation
		while(!seedstack.empty())
		{
			// fetch a seed from the stack
			seed = seedstack.top();
			seedstack.pop();

			// if the seed has been segmented, then go on next
			if (facesegmap.count(seed))
				continue;

			// mark the seed with the segment num
			facesegmap[seed] = regionnum;

			// push it in a new stack
			grow.push(seed);
			// find all of half faces in this segment with this seed
			while(!grow.empty())
			{
				curr = grow.top();
				grow.pop();

				// iterate each half edge of the half face
				for(fe_it = tmesh->face_half_edge_iter(curr); fe_it; ++fe_it)
				{
					// fetch neighbor boundary halfface
					mate = boudary_mate_half_edge_handle(fe_it.handle());
					newface = tmesh->handle_to_entity(mate).half_face_handle();

					// if it has been segmented, do nothing
					if (facesegmap.count(newface))
						continue;

					// if a sharp edge
					if( is_sharpedge(fe_it.handle(), mate))
					{			
						// push it into the seed stack
						seedstack.push(newface);
					}
					// if not a sharp edge
					else
					{
						// segment it and push it into the grow stack
						facesegmap[newface] = regionnum;
						grow.push(newface);
					}
				}
			}
			++ regionnum;
		}
		-- regionnum;
	}


	/* find the mate half edge on the boundary */
	TetraMesh::HalfEdgeHandle VertexSmoother::boudary_mate_half_edge_handle(TetraMesh::HalfEdgeHandle he)
	{
		// check if the index of half edge is valid
		if (he.idx() < 0 || he.idx() > (tmesh->size_halfedge()-1))
		{
			std::cerr << "Half Edge Handle Out of Range!" << std::endl;
			return TetraMesh::HalfEdgeHandle(-1);
		}

		TetraMesh::HalfEdgeHandle mate;

		/* if the tetra contains halfedge has more than one boundary faces, return the mate half face in the same tetra */
		mate = tmesh->mate_half_edge_handle(he);
		if (tmesh->is_boundary(tmesh->handle_to_entity(mate).half_face_handle()))
			return mate;

		/* get one of candidate half edge */
		mate = tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(tmesh->mate_half_edge_handle(he)));
		while (!tmesh->is_boundary(tmesh->handle_to_entity(mate).half_face_handle()) && mate != he && mate.idx() != -1)
		{
			mate = tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(mate));
		}

		if (mate == he)
			mate = TetraMesh::HalfEdgeHandle(-1);
		return mate;
	}

	/* Check if a half edge is sharp
	* he1 and he2 belong to two neighbor halfface respectively, and they share the same edge
	*/
	bool VertexSmoother::is_sharpedge(TetraMesh::HalfEdgeHandle he1, TetraMesh::HalfEdgeHandle he2)
	{
		TetraMesh::HalfFaceHandle hf1, hf2;
		TetraMesh::Normal norm1, norm2;
		double w;

		hf1 = tmesh->handle_to_entity(he1).half_face_handle();
		hf2 = tmesh->handle_to_entity(he2).half_face_handle();

		// fetch face normal
		norm1 = tmesh->normal(hf1);
		norm2 = tmesh->normal(hf2);

		// calculate the normals angle
		w = acos(norm1.normalize() | norm2.normalize());

		if (w > PI/7.2)/* >25бу */
			return true;

		return false;
	}

	/* get vertex information: handle, vertexkind, related vec
	*/
	void VertexSmoother::vertexinformation(SmoothVertex &vinfo, PointHandle p)
	{
		std::vector<TetraHandle> pointStar;
		std::set<int> segments;
		TetraMesh::HedronFaceIter hf_it;
		std::vector<TetraMesh::HalfFaceHandle> boundaryface;

		/* get point handle*/
		vinfo.handle = p;

		// if it is not a boundary point, it is free vertex
		if (!tmesh->is_boundary(p))
		{
			vinfo.kind = FREEVERTEX;
		}
		else
		{
			// fetch tetras incident to the point
			tmesh->point_star(p, pointStar);

			// get the boundary faces incident to the point
			for (unsigned int i = 0; i < pointStar.size(); ++i)
			{
				for (hf_it = tmesh->hedron_face_iter(pointStar[i]); hf_it; ++hf_it)
				{
					if (tmesh->is_boundary(hf_it.handle()) && tmesh->has_point(hf_it.handle(), p))
					{
						// if the halfface is boundary and contains the point
						boundaryface.push_back(hf_it.handle());
						segments.insert(facesegmap[hf_it.handle()]);
					}
				}
			}
			assert(segments.size() > 0);

			// if the point only belongs to one segment, then it is a facet vertex
			if(facetvertex(vinfo))
			{
				vinfo.kind = FACETVERTEX;
			}
			// if the point belongs to two segments, then it is a segment vertex
			else if(segmentvertex(vinfo))
			{
				vinfo.kind = SEGMENTVERTEX;
			}
			// if the point belongs to more than two segments, then it is a fixed vertex
			else
			{
				vinfo.kind = FIXEDVERTEX;
			}
		}
	}

	
	/* optimization-based smoothing for a single vertex 
	 * kinds: smooth type
              SMOOTHFACETVERTICES 0x01
              SMOOTHSEGMENTVERTICES 0x02
              SMOOTHFIXEDVERTICES 0x04
	 */
	bool VertexSmoother::nonsmoothsinglevertex(PointHandle ph, double &worstout, int kinds)
	{
		bool noghosts = true;
		bool smoothed;
		/* a list of all the tets incident to this point */
		std::vector<TetraHandle> pointStar;
		/* a list of information about all the incident tets */
		std::vector<OptTet> incidenttets;
		OptTet opttet;
		Point origpoint;

		/* don't smooth if smoothing disabled */
		if (do_nonsmooth == 0) return false;

		/* get incident tets */
		tmesh->point_star(ph, pointStar);
		assert(pointStar.size() > 0);

		/* copy tags to incident tet data structure */
		for (unsigned int i = 0; i < pointStar.size(); i++)
		{
			/* copy tetrahedron index */
			opttet.handle = pointStar[i];
			incidenttets.push_back(opttet);
		}

		/* smooth the vertex. */
		origpoint = tmesh->point(ph);
		smoothed = nonsmooth(ph, incidenttets, worstout, kinds);
		if (smoothed)
		{
			/* record this vertex insertion in the journal */
			pointvec.clear();thvec.clear();phvec.clear();vhvec.clear();
			pointvec.push_back(origpoint);
			phvec.push_back(ph);
			journals->pushJournal(SMOOTHVERTEX, pointvec, thvec, phvec, vhvec);

			return true;
		}
		return false;
	}


	/* perform non-smooth optimization of a vertex's position.
	*  ph - the point to be smoothed.	 
	*  smoothkinds: smooth type
	                SMOOTHFACETVERTICES 0x01
	                SMOOTHSEGMENTVERTICES 0x02
	                SMOOTHFIXEDVERTICES 0x04
	*  If we end up moving the point, return true.
	*  If for some reason we can't, return false.
	*/
	bool VertexSmoother::nonsmooth(PointHandle ph, std::vector<OptTet> &incidenttets, double &outworstqual, int smoothkinds)
	{
		int numiter = 0;       /* number of optimization iterations */
		double worstqual;      /* the numerical value of the worst quality function */
		double thisqual;       /* quality of the current tet */
		double oldworstqual;   /* the numerical value of the worst quality function at the last step */
		double improvement;    /* how much we've improved this step */
		double dlength;        /* length of d */
		double r;              /* expected rate of improvement */
		double newr;           /* temp r var */
		double alpha;          /* step size */
		double newalpha;       /* candidate new step size */
		double rate;           /* the rate a particular function will improve in the direction d */
		TetraMesh::Normal change;      /* change we'll make to the point */
		Point verts[4];        /* tetra points*/
		SmoothVertex vinfo;    /* contains vertex type and a vector related to the vertex */
		double qe;
		double allgrads[3];
		TetraMesh::Normal d;           /* direction to move vertex */
		/* the gradients in the active set */
		std::vector<TetraMesh::Normal> activegrads;
		Point origpoint;       /* to save the original vertex location */
		Point &p = tmesh->point(ph);   /* the vertex to be altered */
		bool is_nonsmoothed = false;

		/* get point information */
		vertexinformation(vinfo, ph);

		/* check whether to try to smooth facet/segment/fixed vertices */
		if (vinfo.kind == FACETVERTEX && ((smoothkinds & SMOOTHFACETVERTICES) == 0))
			return false;

		if (vinfo.kind == SEGMENTVERTEX && ((smoothkinds & SMOOTHSEGMENTVERTICES) == 0))
			return false;

		if (vinfo.kind == FIXEDVERTEX && ((smoothkinds & SMOOTHFIXEDVERTICES) == 0))
			return false;

		switch (vinfo.kind)
		{
		case FREEVERTEX:
			++ stats->freesmoothattempts;
			break;
		case FACETVERTEX:
			++ stats->facetsmoothattempts;
			break;
		case SEGMENTVERTEX:
			++ stats->segmentsmoothattempts;
			break;
		case FIXEDVERTEX:
			++ stats->fixedsmoothattempts;
			break;
		default:
			printf("I don't know vertex type %d, dying\n", vinfo.kind);
			return false;
		}
		++ stats->nonsmoothattempts;

		/* save the original position of the vertex */
		origpoint = p;

		/* find the worst quality of all incident tets */
		worstqual = HUGEFLOAT; 
		for (unsigned int i = 0; i < incidenttets.size(); i++)
		{
			if (!tmesh->is_valid(incidenttets[i].handle))
				continue;
			tetPoints(tmesh, incidenttets[i].handle, verts);
			thisqual = qualitycal->tetquality(verts[0], verts[1], verts[2], verts[3], qualitymetric);

			/* is this the worst quality we've seen? */
			if (thisqual < worstqual) worstqual = thisqual;
		}
		if (worstqual < 0.0)
			return false;

		/* if we're using surface quadrics */
		if (usequadrics && vinfo.kind == FIXEDVERTEX)
		{
			/* this vertex had better have an initialized quadric */
			assert(quadricVec->quadric(ph).hasquadric);

			/* check whether the surface quadric is the worse than all the tets */
			qe = quadricVec->quadricerrortet(ph);
			if (qe < worstqual) 
				worstqual = qe;
		}

		//outworstqual = worstqual;

		allgrads[0] = allgrads[1] = allgrads[2] = 0.0;

		/* identify the active set A of quality functions that are
		nearly as bad as f, the worst of them all */
		for (unsigned int i = 0; i < incidenttets.size(); i++)
		{
			/* compute gradient info for this tet */
			getoptinfo(&incidenttets[i], vinfo);
		}

		activegrads.clear();
		getactiveset(vinfo, incidenttets, activegrads, worstqual);
		
		if (activegrads.size() == 0)
			return false;

		/* d <- point on the convex hull of all gradients nearest origin */
		minconvexhullpoint(ph, activegrads, d);

		/* if d is the origin, we can't improve this vertex with smoothing. */
		dlength = d.norm();
		if (dlength < DEPSILON)
		{
			return false;
		}

		/* otherwise, it's time to do some smoothing! */
		do
		{
			/* find r, the expected rate of improvement. r is computed as the minimum
			dot product between the improvment direction d and all of the active 
			gradients, so it's like "how fast do we move in the direction of the
			gradient least favored by d." */

			/* start with an absurdly big r */
			r = HUGEFLOAT;
			/* find the smallest dot product */
			for (unsigned int i = 0; i < activegrads.size(); i++)
			{
				newr = d | activegrads[i];

				if (newr <= 0.0)
				{
					/* if we have already moved the vertex some, this is still a success */
					if (p != origpoint)
					{
						outworstqual = worstqual;

						/* record stats */
						switch (vinfo.kind)
						{
						case FREEVERTEX:
							++ stats->freesmoothsuccesses;
							break;
						case FACETVERTEX:
							++ stats->facetsmoothsuccesses;
							break;
						case SEGMENTVERTEX:
							++ stats->segmentsmoothsuccesses;
							break;
						case FIXEDVERTEX:
							++ stats->fixedsmoothsuccesses;
							break;
						default:
							printf("i don't know vertex type %d, dying\n", vinfo.kind);
						}
						++stats->nonsmoothsuccesses;
						return true;
					}
					else
					{
						return false;
					}
				}

				/* this had better be positive */
				assert(newr > 0.0);

				if (newr < r)
				{
					r = newr;
				}
			}

			/* save the worst quality from the previous step */
			oldworstqual = worstqual;

			/* initialize alpha to the nearest intersection with another
			quality function */
			alpha = getinitialalpha(incidenttets, vinfo, d, r, worstqual);

			assert(alpha >= 0.0);

			/* if we didn't find a limit for alpha above, at least limit it
			so that we don't invert any elements. */
			if (alpha == HUGEFLOAT)
			{
				for (unsigned int i = 0; i < incidenttets.size(); i++)
				{
					/* if moving in the direction d will decrease 
					this element's volume */
					rate = d | incidenttets[i].volumegrad;
					if (rate < 0.0)
					{
						newalpha = -incidenttets[i].volume / (2.0 * rate);

						/* if this is smaller than the current step size,
						use it */
						if (newalpha < alpha)
							alpha = newalpha;
					}
				}
			}

			/* do normal line search */
			nonsmoothlinesearch(vinfo, d, worstqual, alpha, r, incidenttets);
			assert(alpha >= 0.0);

			/* move vertex in direction d step size alpha */
			change = TetraMesh::Normal(alpha*d[0], alpha*d[1], alpha*d[2]);
			p = p + change;

			/* recompute quality information */
			oldworstqual = worstqual;
			worstqual = HUGEFLOAT; 
			for (unsigned int i = 0; i < incidenttets.size(); i++)
			{
				if (!tmesh->is_valid(incidenttets[i].handle))
					continue;
				tetPoints(tmesh, incidenttets[i].handle, verts);
				thisqual = qualitycal->tetquality(verts[0], verts[1], verts[2], verts[3], qualitymetric);

				/* is this the worst quality we've seen? */
				if (thisqual < worstqual) worstqual = thisqual;
			}

			/* if we're using surface quadrics */
			if (usequadrics && vinfo.kind == FIXEDVERTEX)
			{
				/* this vertex had better have an initialized quadric */
				assert(quadricVec->quadric(vinfo.handle).hasquadric);

				/* check whether the surface quadric is the worse than all the tets */
				qe = quadricVec->quadricerrortet(ph);
				if (qe < worstqual) 
					worstqual = qe;
			}

			assert(worstqual >= 0.0);

			/* how much did we improve this step? */
			improvement = worstqual - oldworstqual;
			if (improvement > MINSMOOTHITERIMPROVE)
				is_nonsmoothed = true;
				
			/* recompute the active set */
			for (unsigned int i = 0; i < incidenttets.size(); i++)
			{
				if (!tmesh->is_valid(incidenttets[i].handle))
					continue;
				/* compute gradient info for this tet */
				getoptinfo(&incidenttets[i], vinfo);
			}

			activegrads.clear();
			getactiveset(vinfo, incidenttets, activegrads, worstqual);
			assert(activegrads.size() > 0);

			/* d <- point on the convex hull of all gradients nearest origin */
			minconvexhullpoint(ph, activegrads, d);

			numiter++;
			dlength = d.norm();

		} while ((dlength > DEPSILON) && 
			     (numiter < MAXSMOOTHITER) &&
			     (improvement > MINSMOOTHITERIMPROVE));

		/* record this change in the journal */
		//tetPoints(tmesh, incidenttets[0].handle, verts);

		// record the journal
		//insertjournalentry(mesh, SMOOTHVERTEX, verts, 4, origpoint, v);

		outworstqual = worstqual;

		if (!is_nonsmoothed)
			return false;

		/* record stats */
		switch (vinfo.kind)
		{
		case FREEVERTEX:
			++stats->freesmoothsuccesses;
			break;
		case FACETVERTEX:
			++stats->facetsmoothsuccesses;
			break;
		case SEGMENTVERTEX:
			++stats->segmentsmoothsuccesses;
			break;
		case FIXEDVERTEX:
			++stats->fixedsmoothsuccesses;
			break;
		default:
			printf("i don't know vertex type %d, dying\n", vinfo.kind);
		}
		++stats->nonsmoothsuccesses;
		return true;
	}


	/* get the information about this tet needed for non-smooth
	optimization of the current quality measure */
	void VertexSmoother::getoptinfo(OptTet *opttet, SmoothVertex vinfo)
	{
		Point vxtp;                       /* the coordinate of the point */
		Point point[4];                   /* the vertices of the tet */
		double edgelength[3][4];          /* the lengths of each of the edges of the tet */
		TetraMesh::Normal edgegrad[3][4]; /* the gradient of each edge length wrt vtx1 */
		TetraMesh::Normal facenormal[4];  /* the normals of each face of the tet */
		double facearea[4];               /* areas of the faces of the tet */
		TetraMesh::Normal facegrad[4];    /* the gradient of each of the face areas wrt vtx1 */
		double volume;                    /* volume of tetrahedron */
		TetraMesh::Normal volumegrad = TetraMesh::Normal(0.0, 0.0, 0.0);     /* the gradient of the volume of the tet wrt vtx1 */
		int edgecount = 0;                /* keep track of current edge */
		TetraMesh::Normal t, u, v;
		TetraMesh::Normal e1 = TetraMesh::Normal(0.0, 0.0, 0.0);
		TetraMesh::Normal e2 = TetraMesh::Normal(0.0, 0.0, 0.0);

		/* temporary variables */
		Point tempV[4];
		double factor, c, top, bot;
		TetraMesh::Normal diff, term1, term2;
		TetraMesh::Normal gradtop, gradbot, gradquot;

		/* radius ratio vars */
		double Z, twooverZ, gradZtfac, gradZufac, gradZvfac;
		double vdott, tdotu, udotv, tlen2, ulen2, vlen2;
		double faceareasum, rootZareasum;
		double umvlen2, vmtlen2, tmulen2;
		double normfac = sqrt(3.0) * 6.0;
		TetraMesh::Normal gradZt, gradZu, gradZv, gradZ;
		TetraMesh::Normal facegradsum;
		TetraMesh::Normal uminusv, vminust, tminusu;
		TetraMesh::Normal rnrrgradterm1, rnrrgradterm2;
		TetraMesh::Normal E[3] = {TetraMesh::Normal(0.0, 0.0, 0.0), 
			                      TetraMesh::Normal(0.0, 0.0, 0.0), 
					       		  TetraMesh::Normal(0.0, 0.0, 0.0)};

		/* V / lrms^3 ratio vars */
		double edgelengthsum = 0.0;
		double lrms;
		TetraMesh::Normal gradlrms;
		TetraMesh::Normal vlrmsterm1, vlrmsterm2;
		TetraMesh::HalfFaceHandle hf;

		/* get point coordinate */
		vxtp = tmesh->point(vinfo.handle);
		/* get tet vertices */
		tetPoints(tmesh, opttet->handle, tempV);
		/* let tet points set begin with the smoothing point */
		bool is_find = false;
		for (int idx = 0; idx < 4; idx++)
		{
			if (vxtp == tempV[idx])
			{
				point[0] = vxtp;
				point[1] = tempV[(idx+1)%4];
				if (idx & 1)
				{
					point[2] = tempV[(idx+3)%4];
					point[3] = tempV[(idx+2)%4];
				}
				else
				{
					point[2] = tempV[(idx+2)%4];
					point[3] = tempV[(idx+3)%4];
				}
				is_find = true;
				break;
			}
		}
		assert(is_find);

		/* set some vectors */
		t = point[1] - point[0];
		v = point[2] - point[0];
		u = point[3] - point[0];

		/* calculate the volume*6 of the tetrahedron using orientation */
		volume = ((t % v) | u) / 6.0;
		opttet->volume = volume;

		/* for each vertex/face of the tetrahedron */
		for (unsigned int i = 0; i < 4; i++) 
		{
			int j, k, l;
			j = (i + 1) & 3;
			if ((i & 1) == 0) {
				k = (i + 3) & 3;
				l = (i + 2) & 3;
			} else {
				k = (i + 2) & 3;
				l = (i + 3) & 3;
			}
			/* compute the normal for each face */
			/* for each vertex i in the loop, the ith face is the face
			opposite i, so that face's normal is found by taking the
			cross product of two edges of the opposite face */

			/* compute a normal vector to this face */
			facenormal[i] = (point[k] - point[j]) % (point[l] - point[j]);

			/* if i=0, this is also the gradient of the volume * 6
			with respect to vertex 0 */
			if (i == 0)
			{
				opttet->volumegrad = volumegrad = TetraMesh::Normal(-facenormal[i][0]/6.0, -facenormal[i][1]/6.0, -facenormal[i][2]/6.0);
			}

			/* now get the real area */
			opttet->facearea[i] = facearea[i] = sqrt(facenormal[i] | facenormal[i]) / 2.0;

			/* compute the gradient of the area for this face */
			if (i == 0)
			{
				/* this face doesn't include vtx1, gradient is zero */
				opttet->facegrad[i] = facegrad[i] = TetraMesh::Normal(0.0, 0.0, 0.0);
			}
			else
			{
				assert(facearea[i] > 0);
				/* gradient scaled by the face's area */
				factor = 1.0 / (4.0 * facearea[i]);

				/* handle each face separately */
				switch(i)
				{
					/* compute the area of face 1 using u and v */
				case 1:
					e1 = u;
					e2 = v;
					break;
				case 2:
					e1 = t;
					e2 = u;
					break;
				case 3:
					e1 = v;
					e2 = t;
					break;
				}

				/* compute first term of gradient */
				c = e2 | (e1 - e2);
				term1 = TetraMesh::Normal(c*e1[0], c*e1[1], c*e1[2]);

				/* compute the second term */
				c = e1 | (e1 - e2);
				term2 = TetraMesh::Normal(c*e2[0], c*e2[1], c*e2[2]);

				/* now, combine the terms, scaled with the 1/4A */
				diff = term1 - term2;
				opttet->facegrad[i] = facegrad[i] = TetraMesh::Normal(factor*diff[0], factor*diff[1], factor*diff[2]);
			}


			/* compute edge lengths for quality measures that need them */
			if (qualitymetric == QUAL_MINSINE ||
				qualitymetric == QUAL_WARPEDMINSINE ||
				qualitymetric == QUAL_VLRMS3RATIO)
			{
				for (int j = i + 1; j < 4; j++) 
				{
					/* e1 is edge from point i to point j */
					e1 = point[j] - point[i];
					opttet->edgelength[i][j] = edgelength[i][j] = e1.norm();

					/* also compute the gradient of the length of this edge */

					/* if vtx1 isn't one of the edge's endpoints, the gradent is zero */
					if (i != 0)
					{
						opttet->edgegrad[i][j] = edgegrad[i][j] = TetraMesh::Normal(0.0, 0.0, 0.0);
					}
					/* otherwise, it's in the negative direction of this edge,
					and scaled by edge length */
					else
					{ 
						factor = -1.0 / edgelength[i][j];
						edgegrad[i][j] = TetraMesh::Normal(factor*e1[0], factor*e1[1], factor*e1[2]);
						opttet->edgegrad[i][j] = edgegrad[i][j];
					}
				}
			}
		}

		/* if the quality measure is minimum sine */
		if ((qualitymetric == QUAL_MINSINE) || (qualitymetric == QUAL_WARPEDMINSINE))
		{
			/* for each edge in the tetrahedron */
			int i, j, k, l;
			for (i = 0; i < 3; i++) 
			{
				for (j = i + 1; j < 4; j++)
				{
					k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
					l = 6 - i - j - k;

					/* compute the sine of this dihedral angle */
					opttet->sine[edgecount] = (3 * volume * edgelength[i][j]) / (2 * facearea[k] * facearea[l]);

					/* if we are warping the minimum sine */
					if (qualitymetric == QUAL_WARPEDMINSINE)
					{
						/* and this is an obtuse angle */
						if ((facenormal[k] | facenormal[l]) > 0)
						{
							/* scale the sin down by WARPFACTOR */
							opttet->sine[edgecount] *= sinewarpfactor;
						}
					}

					/* compute the gradient of the sine           3 * V * lij
					we need the gradient of this expression:     -------------
					                                              2 * Ak * Al
					so, find the gradient of the top product, the bottom product, then the quotient
					*/
					top = volume * edgelength[i][j];
					bot = facearea[k] * facearea[l];

					/* find gradient of top */
					gradtop = gradproduct(volume, edgelength[i][j], volumegrad, edgegrad[i][j]);

					/* find gradient of bottom */
					gradbot = gradproduct(facearea[k], facearea[l], facegrad[k], facegrad[l]);

					/* now, find the gradient of the quotient */
					gradquot = gradquotient(top, bot, gradtop, gradbot);

					/* now scale with constant factor */
					opttet->sinegrad[edgecount] = TetraMesh::Normal(1.5*gradquot[0], 1.5*gradquot[1], 1.5*gradquot[2]);

					/* if this is a facet vertex, project gradient onto facet */
					if (vinfo.kind == FACETVERTEX)
					{
						opttet->sinegrad[edgecount] = vprojecttoplane(opttet->sinegrad[edgecount], vinfo.vec);
					}

					/* if this is a segment vertex, project gradient onto segment */
					if (vinfo.kind == SEGMENTVERTEX)
					{
						opttet->sinegrad[edgecount] = vproject(opttet->sinegrad[edgecount], vinfo.vec);
					}

					edgecount++;
				}
			}
		}

		/* compute stuff for radius ratio */
		if (qualitymetric == QUAL_RADIUSRATIO)
		{
			/* compute intermediate quantity Z */
			Z = getZ(point[0], point[1], point[2], point[3]);
			twooverZ = 2.0 / Z;

			/* some dot products */
			vdott = v | t;
			tdotu = t | u;
			udotv = u | v;

			/* some vector lengths */
			uminusv = u - v;
			vminust = v - t;
			tminusu = t - u;
			tlen2 = (t[0] * t[0]) + (t[1] * t[1]) + (t[2] * t[2]);
			ulen2 = (u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]);
			vlen2 = (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]);
			umvlen2 = (uminusv[0] * uminusv[0]) + (uminusv[1] * uminusv[1]) + (uminusv[2] * uminusv[2]);
			vmtlen2 = (vminust[0] * vminust[0]) + (vminust[1] * vminust[1]) + (vminust[2] * vminust[2]);
			tmulen2 = (tminusu[0] * tminusu[0]) + (tminusu[1] * tminusu[1]) + (tminusu[2] * tminusu[2]);

			/* compute Z's gradient */
			gradZtfac = twooverZ * 
				(
				(ulen2 * vdott - vlen2 * tdotu) * (ulen2 - vlen2) - 
				(ulen2 * vlen2 + tlen2 * udotv) * (umvlen2)
				);
			gradZufac = twooverZ * 
				(
				(vlen2 * tdotu - tlen2 * udotv) * (vlen2 - tlen2) - 
				(vlen2 * tlen2 + ulen2 * vdott) * (vmtlen2)
				);
			gradZvfac = twooverZ * 
				(
				(tlen2 * udotv - ulen2 * vdott) * (tlen2 - ulen2) - 
				(tlen2 * ulen2 + vlen2 * tdotu) * (tmulen2)
				);

			/* compute t, u, v components of gradient */
			gradZt = TetraMesh::Normal(gradZtfac*t[0], gradZtfac*t[1], gradZtfac*t[2]);
			gradZu = TetraMesh::Normal(gradZufac*u[0], gradZufac*u[1], gradZufac*u[2]);
			gradZv = TetraMesh::Normal(gradZvfac*v[0], gradZvfac*v[1], gradZvfac*v[2]);

			/* add the components together to form grad(Z) */
			gradZ = gradZt + gradZu + gradZv;

			/* compute sqrt (Z * (sum of face areas)) */
			faceareasum = opttet->facearea[0] +
				          opttet->facearea[1] +
				          opttet->facearea[2] +
				          opttet->facearea[3];
			rootZareasum = sqrt(Z * faceareasum);

			/* set the actual root normalized radius ratio */
			opttet->rnrr = (normfac * volume) / rootZareasum;

			assert(opttet->rnrr > 0.0);

			/* sum of face gradients */
			facegradsum = facegrad[0] + facegrad[1] + facegrad[2] + facegrad[3];

			/* compute the first term */
			rnrrgradterm1 = TetraMesh::Normal((1.0/rootZareasum)*volumegrad[0], 
				                              (1.0/rootZareasum)*volumegrad[1],
											  (1.0/rootZareasum)*volumegrad[2]);

			/* compute the second term */
			facegradsum = TetraMesh::Normal(Z*facegradsum[0], Z*facegradsum[1], Z*facegradsum[2]);
			gradZ = TetraMesh::Normal(faceareasum*gradZ[0], faceareasum*gradZ[1], faceareasum*gradZ[2]);
			rnrrgradterm2 = facegradsum + gradZ;
			c = volume / (2 * (rootZareasum * rootZareasum * rootZareasum));
			rnrrgradterm2 = TetraMesh::Normal(c*rnrrgradterm2[0], c*rnrrgradterm2[1], c*rnrrgradterm2[2]);

			/* finally, compute the gradient of the radius ratio */
			opttet->rnrrgrad = rnrrgradterm1 - rnrrgradterm2;
			opttet->rnrrgrad = TetraMesh::Normal(normfac*opttet->rnrrgrad[0], normfac*opttet->rnrrgrad[1], normfac*opttet->rnrrgrad[2]);

			/* if this is a facet vertex, project gradient onto facet */
			/*if (vtxkind == FACETVERTEX || vtxkind == FIXEDVERTEX)*/
			if (vinfo.kind == FACETVERTEX)
			{
				opttet->rnrrgrad = vprojecttoplane(opttet->rnrrgrad, vinfo.vec);
			}

			/* if this is a segment vertex, project gradient onto segment */
			if (vinfo.kind == SEGMENTVERTEX)
			{
				opttet->rnrrgrad = vproject(opttet->rnrrgrad, vinfo.vec);
			}
		}

		/* if the quality measure is volume to edge length ratio */
		if (qualitymetric == QUAL_VLRMS3RATIO)
		{
			/* for each edge in the tetrahedron */
			int i, j, k, l;
			for (i = 0; i < 3; i++) 
			{
				for (j = i + 1; j < 4; j++) 
				{
					k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
					l = 6 - i - j - k;

					/* accumulate edge length sum */
					edgelengthsum += edgelength[i][j] * edgelength[i][j];
				}
			}

			/* compute the root mean square */
			lrms = sqrt((1.0 / 6.0) * edgelengthsum);

			normfac = 6.0 * sqrt(2.0);

			/* compute the raw ratio */
			opttet->vlrms3r = (normfac * volume) / (lrms * lrms * lrms);

			/* compute gradient of lrms */
			gradlrms = t + u + v;
			c = (-1.0 / (6.0 * lrms));
			gradlrms = TetraMesh::Normal(c*gradlrms[0], c*gradlrms[1], c*gradlrms[2]);

			/* compute the terms of the gradient of the ratio */
			c = 1.0 / (lrms * lrms * lrms);
			vlrmsterm1 = TetraMesh::Normal(c*volumegrad[0], c*volumegrad[1], c*volumegrad[2]);
			c = (3.0 * volume) / (lrms * lrms * lrms * lrms);
			vlrmsterm2 = TetraMesh::Normal(c*gradlrms[0], c*gradlrms[1], c*gradlrms[2]);

			/* add terms and normalize */
			opttet->vlrms3rgrad = vlrmsterm1 - vlrmsterm2;
			opttet->vlrms3rgrad = TetraMesh::Normal(normfac*opttet->vlrms3rgrad[0],
				                                    normfac*opttet->vlrms3rgrad[1],
													normfac*opttet->vlrms3rgrad[2]);

			/* if this is a facet or fixed vertex, 
			project gradient onto facet or plane */
			/*if (vtxkind == FACETVERTEX || vtxkind == FIXEDVERTEX)*/
			if (vinfo.kind == FACETVERTEX)
			{
				opttet->vlrms3rgrad = vprojecttoplane(opttet->vlrms3rgrad, vinfo.vec);
			}

			/* if this is a segment vertex, project gradient onto segment */
			if (vinfo.kind == SEGMENTVERTEX)
			{
				opttet->vlrms3rgrad = vproject(opttet->vlrms3rgrad, vinfo.vec);
			}
		}
	}


	/* compute Z, a quantity associated with circumradius computation
	TODO this code is lifted from Jonathan's tetcircumcenter computation
	in primitives.c */
	double VertexSmoother::getZ(TetraMesh::Normal tetorg,
		                        TetraMesh::Normal tetdest,
		                        TetraMesh::Normal tetfapex,
		                        TetraMesh::Normal tettapex)
	{
		TetraMesh::Normal ot, dt, ft;
		double otlength, dtlength, ftlength;
		TetraMesh::Normal crossdf, crossfo, crossod;
		TetraMesh::Normal ct;

		/* Use coordinates relative to the apex of the tetrahedron. */
		ot = tetorg - tettapex;
		dt = tetdest - tettapex;
		ft = tetfapex - tettapex;
		/* Squares of lengths of the origin, destination, and face apex edges. */
		otlength = ot | ot;
		dtlength = dt | dt;
		ftlength = ft | ft;
		/* Cross products of the origin, destination, and face apex vectors. */
		crossdf = dt % ft;
		crossfo = ft % ot;
		crossod = ot % dt;

		/* Calculate offset (from apex) of circumcenter. */
		ct[0] = (otlength * crossdf[0] + dtlength * crossfo[0] + ftlength * crossod[0]);
		ct[1] = (otlength * crossdf[1] + dtlength * crossfo[1] + ftlength * crossod[1]);
		ct[2] = (otlength * crossdf[2] + dtlength * crossfo[2] + ftlength * crossod[2]);

		/* Calculate the length of this vector, which is Z */
		return sqrt(ct | ct);
	}


	/* given a set of tets incident to a vertex, and the quality
	of the worst quality function that varies with that vertex,
	compute the active set A of quality functions very near
	the worst */

	void VertexSmoother::getactiveset(SmoothVertex sv, std::vector<OptTet> &incidenttets, 
		                              std::vector<TetraMesh::Normal> &activegrads, double worstqual)
	{
		/* if we are including surface quadrics, give them a chance to
		enter the active set */
		assert(incidenttets.size() >= 0);

		/* if (improvebehave.usequadrics && hasquadric(vtx) && ((struct vertextype *) arraypoolforcelookup(&vertexinfo, vtx))->kind != FREEVERTEX) */
		if (usequadrics && quadricVec->quadric(sv.handle).hasquadric && sv.kind == FIXEDVERTEX)
		{
			/* if the quadric is close enough to worst, add it to the active set */
			if (quadricVec->quadricerrortet(sv.handle) <= (worstqual * ACTIVESETFACTOR))
			{
				/* fetch the gradient of this quadric */
				TetraMesh::Normal quadgrad;
				quadgrad = quadricVec->quadricgradtet(sv.handle);

				/* copy the gradient into the list of active gradients */
				activegrads.push_back(quadgrad);
			}
		}

		for (unsigned int i = 0; i < incidenttets.size(); i++)
		{
			if (!tmesh->is_valid(incidenttets[i].handle))
				continue;

			switch (qualitymetric)
			{
			case QUAL_WARPEDMINSINE:
			case QUAL_MINSINE:
				for (int j = 0; j < 6; j++)
				{
					/* is this close enough to the worst? */
					if (incidenttets[i].sine[j] < (worstqual * ACTIVESETFACTOR))
					{
						/* get the actual gradient value */
						activegrads.push_back(incidenttets[i].sinegrad[j]);
					}
				}
				break;
			case QUAL_RADIUSRATIO:

				if (incidenttets[i].rnrr <= (worstqual * ACTIVESETFACTOR))
				{
					/* get the actual gradient value */
					activegrads.push_back(incidenttets[i].rnrrgrad);
				}
				break;
			case QUAL_VLRMS3RATIO:

				if (incidenttets[i].vlrms3r <= (worstqual * ACTIVESETFACTOR))
				{
					/* get the actual gradient value */
					activegrads.push_back(incidenttets[i].vlrms3rgrad);
				}
				break;
			default:
				printf("i don't know how to compute the active set for qual measure %d\n", qualitymetric);
				break;
			}
		}

		/* we must have at least found the worst */
		/* assert(*numactive > 0); */
	}



	/* finds the point on the convex hull of P nearest the origin */
	void VertexSmoother::minconvexhullpoint(PointHandle ph, std::vector<TetraMesh::Normal> P, TetraMesh::Normal &nearest)
	{
		std::vector<TetraMesh::Normal> B;         /* the basis for the convex hull point */
		std::vector<TetraMesh::Normal> empty;     /* empty set for known basis */
		Point p, q, r;                            /* basis points */
		TetraMesh::Normal pmq;                    /* p minus q */
		double c, d, l;                           /* scalar factors */
		TetraMesh::Normal s, t, s2, t2;           /* temporary points */
		TetraMesh::Normal sxt, sxr, rxt;          /* temporary cross products */

		assert(P.size() > 0);

		/* find a basis for the minimum point on the convex hull */
		findbasis(P, empty, B);

		switch(B.size())
		{
			/* if the basis is just a single point, return that point */
		case 1:
			nearest = B[0];
			return;
			break;
			/* for two points, find the closest point to the origin on
			the line between the two points */
		case 2:
			p = B[0];
			q = B[1];
			pmq = p - q;

			/*
			nearest = q - dot(q,p-q)/(length(p-q)^2) * (p-q)
			*/
			l = pmq.norm();

			/* if these points are the same, just return one of them */
			if (l == 0.0)
			{
				nearest = B[0];
				return;
			}

			c = (q | pmq) / (l * l);
			nearest = TetraMesh::Normal(c*pmq[0], c*pmq[1], c*pmq[2]);
			nearest = q - nearest;

			return; 
			break;
			/* for three points, find the point closest to the origin
			on the triangle that the they form */
		case 3:
			p = B[0];
			q = B[1];
			r = B[2];

			s = p - r;
			t = q - r;

			sxt = s % t;
			rxt = r % t;
			sxr = s % r;

			/* if any of these cross products is really tiny, give up
			and return the origin */
			if (sxt.norm() < NEARESTMIN || rxt.norm() < NEARESTMIN || sxr.norm() < NEARESTMIN)
			{
				nearest = TetraMesh::Normal(0.0, 0.0, 0.0);
				return;
			}

			c = (sxt | rxt) / (sxt | sxt);
			d = (sxt | sxr) / (sxt | sxt);

			s2 = TetraMesh::Normal(c*s[0], c*s[1], c*s[2]);
			t2 = TetraMesh::Normal(d*t[0], d*t[1], d*t[2]);

			nearest = r - s2 - t2;

			return;
			break;
			/* if the basis has four points, they must enclose the origin
			so just return the origin. */
		case 4:
			nearest = TetraMesh::Normal(0.0, 0.0, 0.0);
			return;
			break;
		default:
			printf("Error, basis size %d is bogus, dying\n", B.size());
			break;
		}
	}



	/* returns the basis B of S union M. S is a set of points known
	to be in the basis */
	void VertexSmoother::findbasis(std::vector<TetraMesh::Normal> M,  std::vector<Point> &S, std::vector<Point> &B)
	{
		Point p, q, r, s;                          /* vertices */
		Point s1, t1, d1, d2;                      /* temporary vertices */
		Point origin = Point(0.0, 0.0, 0.0);
		std::vector<Point> localS;                 /* for passing to recursive calls */
		std::vector<Point> localM;                 /* for passing to recursive colls */

		assert(M.size() > 0);

		/* we assume that M was passed to us shuffled, so that taking
		the last element is equivalent to removing a random one. */
		p = M[M.size()-1];

		/* if M has only one element */
		if (M.size() == 1)
		{
			/* and S has no elements */
			if (S.size() == 0)
			{
				/* then the single element in M must be the
				entire basis, just send back M */
				B.push_back(M[0]);
				return;
			}

			/* otherwise, because we assume the last element
			we just removed is not part of the basis, assign
			the basis to be the elements of S */
			copy(S.begin(), S.end(), back_inserter(B));
		}
		/* M has more than one element. Throw one out (p), and look for the
		basis assuming p is not part of it. */
		else
		{
			/* make a new copy of M minus the last element */
			copy(M.begin(), M.end()-1, back_inserter(localM));
			copy(S.begin(), S.end(), back_inserter(localS));
			findbasis(localM, localS, B);
		}

		/* now the we have determined the basis without p, we need to 
		go back and check whether p actually is part of the basis. */

		switch (B.size())
		{
			/* if the returned basis has just one point q, we just need to check
			whether p is closer to the origin than q */
		case 1:
			/* fetch the actual coordinates from the mesh */
			q = B[0];

			/* compute the vector from q to p */
			d1 = p - q;

			/* check the sign of the dot product. >=0 means p doesn't
			improve the basis */
			if ( (q | d1) >= 0)
			{
				/* this is a good B, send it back!*/
				return;
			}
			break;

			/* check whether p improves the basis using math I don't understand */        
		case 2:
			/* fetch coordinates from the mesh */
			q = B[0];
			r = B[1];

			/* compute vector s from r to p */
			s1 = p - r;
			/* compute vector t from r to q */
			t1 = q - r;

			/* now a couple of cross products */
			d1 = s1 % t1;
			d2 = r % t1;

			/* can p improve the basis? */
			if ((d1 | d2) >= 0)
			{
				/* nope! send back B as is. */
				return;
			}            
			break;
		case 3:
			/* fetch coordinates from the mesh */
			q = B[0];
			r = B[1];
			s = B[2];

			/* does p improve the basis? */
			if (orient(p, q, r, s) * orient(origin, q, r, s) <= 0)
			{
				/* nope! send back B as is. */
				return;
			}
			break;
		default:
			/* B has size of 4, and any basis of this size is optimal */
			return;
			break;
		}

		/* if we have made it this far, we know that p actually is a part of
		any basis of S union M */

		/* if p was the last element of M, or if S already has three other basis
		points, we're done and just need to send back S union p. */
		if ((M.size() == 1) || (S.size() == 3))
		{
			/* copy S into B */
			B.clear();
			copy(S.begin(), S.end(), back_inserter(B));
			/* and add p at the end */
			B.push_back(p);

			return;
		}
		/* there may be more basis points to find! move p from M to the known
		basis point set S, and go again */
		else
		{
			/* create the new S */
			localS.clear();
			copy(S.begin(), S.end(), back_inserter(localS));
			/* add p to the end of it */
			localS.push_back(p);

			/* create the new M, leaving off the last element */
			localM.clear();
			copy(M.begin(), M.end()-1, back_inserter(localM));

			/* find any basis points remaining in M */
			findbasis(localM, localS, B);

			return;
		}        
	}

	/* for our initial step size, we use the distance
	to the next intersection with another quality funtion.
	this is the point at which the other quality function
	becomes the worst. we use a single-term taylor expansion
	to approximate all of the quality functions as lines,
	so we'll have to do a line search to find our ultimate
	step size. */
	double VertexSmoother::getinitialalpha(std::vector<OptTet> &incidenttets, SmoothVertex vinfo, 
		                                   TetraMesh::Normal d, double r, double worstqual)
	{
		double alpha = HUGEFLOAT;
		double newalpha;
		double rate;

		/* if we are including surface quadrics add check for it */
		assert(incidenttets.size() > 0);
		if (usequadrics && quadricVec->quadric(vinfo.handle).hasquadric && vinfo.kind == FIXEDVERTEX)
		{
			/* fetch the gradient of this quadric */
			TetraMesh::Normal quadgrad;
			quadgrad = quadricVec->quadricgradtet(vinfo.handle);

			/* if this function improves more slowly than any in the active set, 
			   then it might end up as the objective. */
			rate = d | quadgrad;
			if (rate + RATEEPSILON < r)
			{
				/* compute the approximation of when this function will become the objective */
				newalpha = (quadricVec->quadricerrortet(vinfo.handle) - worstqual) / (r - rate);

				/* if this is smaller than our current step size, use it for the step size */
				if (newalpha < alpha)
					alpha = newalpha;
			}
		}

		for (unsigned int i = 0; i < incidenttets.size(); i++)
		{
			switch (qualitymetric)
			{
			case QUAL_WARPEDMINSINE:
			case QUAL_MINSINE:
				for (int j = 0; j < 6; j++)
				{
					/* if this function improves more slowly than any in the active set, 
					   then it might end up as the objective. */
					rate = d | incidenttets[i].sinegrad[j];
					if (rate + RATEEPSILON < r)
					{
						/* compute the approximation of when this
						   function will become the objective */
						newalpha = (incidenttets[i].sine[j] - worstqual) / (r - rate);

						/* if this is smaller than our current step size, use it for the step size */
						if (newalpha < alpha)
							alpha = newalpha;
					}
				}
				break;
			case QUAL_RADIUSRATIO:
				/* if this function improves more slowly than any in the active set, 
				   then it might end up as the objective. */
				rate = d | incidenttets[i].rnrrgrad;
				if (rate + RATEEPSILON < r)
				{
					/* compute the approximation of when this function will become the objective */
					newalpha = (incidenttets[i].rnrr - worstqual) / (r - rate);

					/* if this is smaller than our current step size, use it for the step size */
					if (newalpha < alpha)
					{
						alpha = newalpha;
					}
				}
				break;
			case QUAL_VLRMS3RATIO:
				/* if this function improves more slowly
				than any in the active set, then it might
				end up as the objective. */
				rate = d | incidenttets[i].vlrms3rgrad;
				if (rate + RATEEPSILON < r)
				{
					/* compute the approximation of when this
					function will become the objective */
					newalpha = (incidenttets[i].vlrms3r - worstqual) / (r - rate);

					/* if this is smaller than our current step size, use it for the step size */
					if (newalpha < alpha)
						alpha = newalpha;
				}
				break;
			default:
				printf("i don't know how to compute alpha for qual measure %d\n", qualitymetric);
				break;
			}
		}
		if (alpha < 0.0)
		{
			alpha = 0.0;
		}
		return alpha;
	}


	/* find the best step to take to improve all of the quality
	functions affected by a vertex vtx in the search direction
	d with an expected rate of improvement r */
	void VertexSmoother::nonsmoothlinesearch(SmoothVertex vinfo, TetraMesh::Normal d, double inworstqual, double & alpha,
		                                     double r, std::vector<OptTet> &incidenttets)
	{
		int numiter = 0;              /* number of iterations */
		Point *v;                     /* the vertex to be modified */
		Point origvertex;             /* to save the original vertex position */
		TetraMesh::Normal offset;     /* the offset to move the vertex, alpha * d */
		double worstqual;             /* the current worst quality */
		double origworstqual;         /* the original worst quality */
		double thisqual;              /* the quality of the current tet */
		double oldworstqual;          /* the worst quality at the last step */
		Point tetv[4];

		/* save the original worst quality */
		origworstqual = oldworstqual = inworstqual;

		/* fetch the original vertex coordinates from the mesh */
		v = &tmesh->point(vinfo.handle);
		origvertex = *v;

		/* keep trying until alpha gets too small or we exceed maximum
		number of iterations */
		while ((alpha > MINSTEPSIZE) && (numiter < MAXLINEITER))
		{
			/* compute the offset from original vertex positon,
			alpha * d */
			offset = TetraMesh::Normal(alpha*d[0], alpha*d[1], alpha*d[2]);
			/* move the vertex */
			*v = *v + offset;

			/* recompute all of the quality functions affected
			by v's position, taking note of the smallest one */
			worstqual = HUGEFLOAT; 
			for (unsigned int i = 0; i < incidenttets.size(); i++)
			{
				tetPoints(tmesh, incidenttets[i].handle, tetv);
				thisqual = qualitycal->tetquality(tetv[0], tetv[1], tetv[2], tetv[3], qualitymetric);

				/* is this the worst quality we've seen? */
				if (thisqual < worstqual) worstqual = thisqual;
			}

			/* if we're using surface quadrics */
			if (usequadrics && quadricVec->quadric(vinfo.handle).hasquadric && vinfo.kind == FIXEDVERTEX)
			{
				/* check whether the surface quadric is the worse than all the tets */
				double qe = quadricVec->quadricerrortet(vinfo.handle);
				if (qe < worstqual) 
					worstqual = qe;
			}

			/* if this is not the first iteration, and
			we did better on the last iteration, use
			the step size from the previous iteration */
			if ((oldworstqual > origworstqual) && (oldworstqual > worstqual))
			{
				/* use the previous step's alpha */
				alpha = (alpha) * 2;
				assert(alpha > 0.0);

				/* put vertex back where it started */
				*v = origvertex;
				return;
			}

			/* if we have succeeded in gaining 90% of the expected
			improvement, accept this initial step size */
			if ((worstqual - origworstqual) > (0.9 * (alpha) * r))
			{
				/* put vertex back where it started */
				*v = origvertex;
				return;
			}

			/* put vertex back where it started */
			*v = origvertex;

			/* cut alpha down by half and try again */
			alpha = (alpha) / 2.0;

			/* save the worst quality from this step */
			oldworstqual = worstqual;
		}

		/* no positive alpha could be found that improved things... give up and return zero */
		alpha = 0.0;
	}

	/* given a boundary vertex and a list of all of it's
	incident tets, determine if it lies in an input facet */
	bool VertexSmoother::facetvertex(SmoothVertex &vinfo)
	{
		int numfaces=0;                         /* total coplanar faces surrounding vertex */
		TetraMesh::Normal refnorm;              /* "reference" normal of first face */
		TetraMesh::Normal curnorm;              /* the current face's normal */
		double dotprod=0.0;                     /* dot product of current and reference normal */
		std::vector<TetraHandle> pointStar;     /* list of tets incident to the point */
		std::vector<TetraMesh::HalfFaceHandle> boundfacelist; /* list of boundary faces for a particular tet */
		TetraMesh::HedronFaceIter hf_iter;

		/* first, check to see if we lie in a facet. to do this, we need
		to build a list of all the boundary faces incident to this vertex.
		this is the same as getting all of the boundary faces of all incident
		tets that include this vertex. */

		/* get incident tet */
		tmesh->point_star(vinfo.handle, pointStar);

		/* get all the boundary faces contains the vertex */
		for (unsigned int i = 0; i < pointStar.size(); ++i)
		{
			for (hf_iter = tmesh->hedron_face_iter(pointStar[i]); hf_iter; ++hf_iter)
			{
				if (tmesh->is_boundary(hf_iter.handle()) && tmesh->has_point(hf_iter.handle(), vinfo.handle))
					boundfacelist.push_back(hf_iter.handle());
			}
		}

		/* for each boundary face*/
		for (unsigned int i = 0; i < boundfacelist.size(); ++i)
		{
			++ numfaces;

			/* find the (unit-length) face normal of this face */
			curnorm = tmesh->update_face_normal(boundfacelist[i]);
			curnorm.normalize();

			/* if this is the first face we have found, establish it's
			normal as the "reference normal" which future faces
			will be tested against for coplanarity */
			if (i == 0)
			{
				refnorm = curnorm;
			}
			else
			{
				/* is the dot product of this normal within
				COPLANARTOL of the reference normal? */
				dotprod = curnorm | refnorm;

				if ((dotprod > (1.0 + COPLANARTOL)) || (dotprod < (1.0 - COPLANARTOL)))
				{
					/* this new normal is not close enough to coplanar.
					we know that this vertex cannot be a facet vertex */
					return false;
				}
			}
		}

		/* if we get here having seen at least 3 coplanar faces, must be facet vertex */
		if (numfaces >= 3)
		{
			vinfo.vec = refnorm;
			return true;
		}
		return false;
	}

	/* given a boundary vertex and a list of all of it's
	incident tets, determine if it lies in an input segment */
	bool VertexSmoother::segmentvertex(SmoothVertex &vinfo)
	{
		bool foundface=false;                   /* have we encountered an incident boundary face yet? */
		TetraMesh::Normal refnorm1;             /* "reference" normal of first face */
		TetraMesh::Normal refnorm2;             /* the other reference plane */
		TetraMesh::Normal curnorm;              /* the current face's normal */
		double dotprod=0.0;                     /* dot product of current and reference normal */
		std::vector<TetraHandle> pointStar;     /* list of tets incident to the point */
		std::vector<TetraMesh::HalfFaceHandle> boundfacelist; /* list of boundary faces for a particular tet */
		int numfaces=0;                                       /* total coplanar faces surrounding vertex */
		TetraMesh::HedronFaceIter hf_iter;

		/* first, check to see if we lie in a facet. to do this, we need
		to build a list of all the boundary faces incident to this vertex.
		this is the same as getting all of the boundary faces of all incident
		tets that include this vertex. */

		/* get incident tet */
		tmesh->point_star(vinfo.handle, pointStar);

		/* get all the boundary faces contains the vertex */
		for (unsigned int i = 0; i < pointStar.size(); ++i)
		{
			for (hf_iter = tmesh->hedron_face_iter(pointStar[i]); hf_iter; ++hf_iter)
			{
				if (tmesh->is_boundary(hf_iter.handle()) && tmesh->has_point(hf_iter.handle(), vinfo.handle))
					boundfacelist.push_back(hf_iter.handle());
			}
		}


		/* for each boundary face*/
		for (unsigned int i = 0; i < boundfacelist.size(); ++i)
		{
			++ numfaces;

			/* find the (unit-length) face normal of this face */
			curnorm = tmesh->update_face_normal(boundfacelist[i]);
			curnorm.normalize();

			/* if this is the first face we have found, establish it's
			normal as the "reference normal" which future faces
			will be tested against for coplanarity */
			if (i == 0)
			{
				refnorm1 = curnorm;
			}
			else
			{
				/* is the dot product of this normal within
				COPLANARTOL of the first reference normal? */
				dotprod = curnorm | refnorm1;

				/* do dot products match for first normal? */
				if ((dotprod > (1.0 + COPLANARTOL)) || (dotprod < (1.0 - COPLANARTOL)))
				{
					/* nope. is there a second reference yet? */
					if (!foundface)
					{
						/* no second reference yet. make this it */
						foundface = true;
						refnorm2 = curnorm;
					}
					else
					{
						dotprod = curnorm | refnorm2;
						if ((dotprod > (1.0 + COPLANARTOL)) || (dotprod < (1.0 - COPLANARTOL)))
							/* this normal doesn't match either of the references */
							return false;
					}
				}
			}
		}

		/* if we get here having seen at least 3 faces, must be segment vertex */
		if (numfaces >= 3)
		{
			/* the segment that this vertex lies on is perpendicular to both reference normals */
			vinfo.vec = refnorm1 % refnorm2;
			return true;
		}
		return false;
	}


	/* perform a pass of combination Laplacian / optimization-based smoothing. */
	bool VertexSmoother::smoothpass(std::vector<PointHandle> *points, double threshold, 
		                            double bestmeans[], double meanqualafter[], 
									double *minqualafter, int smoothkinds)
	{
		/* check mesh */
		if (tmesh == NULL)
			return false;
		/* initialize */
		initialize();

		int optattempts=0;             /* number of times optimization-based is tried */
		int optsuccesses=0;            /* number of times it succeeds */
		int lapattempts=0;
		int lapsuccesses=0;
		int lapverts[2] = {0,0};
		int fixedverts[2] = {0,0};
		int facetverts[2] = {0,0};
		int segmentverts[2] = {0,0};
		int freeverts[2] = {0,0};
		double worstqual;
		bool smoothed;                 /* was smoothing successful? */
		int kind;                      /* number of degrees of freedom for vertex */
		double minqualbefore = HUGEFLOAT;
		double origqual;
		//double meanqualbefore[NUMMEANTHRESHOLDS];
		int nonexist = 0;
		SmoothVertex smoothedvert;
		std::vector<SmoothVertex> smoothedverts; /* list of vertices that have already been smoothed */
		bool skip = false;
		bool dynfailcondition = true;
		std::set<TetraHandle> tetraset;
		std::vector<TetraHandle> tetravec;

		if (points != NULL)
		{
			/* get all the tetras incident to the smoothing points */
			for(unsigned int pointidx = 0; pointidx < points->size(); ++pointidx)
			{
				tetravec.clear();
				tmesh->point_star((*points)[pointidx], tetravec);
				for (unsigned int i = 0; i < tetravec.size(); i++)
					tetraset.insert(tetravec[i]);
			}
			tetravec.clear();
			copy(tetraset.begin(), tetraset.end(), back_inserter(tetravec));

			/* calculate the original quality */
			origqual = qualitycal->minstackquality(tmesh, tetravec, qualitymetric);
		}
		else
			origqual = qualitycal->minmeshquality(tmesh, qualitymetric);

		/* do vertex smoothing */
		do 
		{
			/* fetch the worst quality before a new pass of smoothing */
			if (points != NULL)
				minqualbefore = qualitycal->minstackquality(tmesh, tetravec, qualitymetric);
			else
				minqualbefore = qualitycal->minmeshquality(tmesh, qualitymetric);

			*minqualafter = HUGEFLOAT;

			/* try to smooth all the points int the list */
			int pointcnt;
			PointHandle currtpoint;
			if (points != NULL)
				pointcnt = points->size();
			else
				pointcnt = tmesh->size_point();
			for (unsigned int pointidx = 0; pointidx < pointcnt; ++pointidx)
			{
				if (points != NULL)
					currtpoint = (*points)[pointidx];
				else
					currtpoint = PointHandle(pointidx);

				/* make sure this point exists */
				if (!tmesh->is_valid(currtpoint))
				{
					nonexist++;
					continue;
				}

				/* look up vertex type */
				vertexinformation(smoothedvert, currtpoint);
				kind = smoothedvert.kind;
				/* record smoothing attempt */
				switch (kind)
				{
				case FREEVERTEX: 
					freeverts[0]++;
					break;
				case FACETVERTEX:
					facetverts[0]++;
					break;
				case SEGMENTVERTEX:
					segmentverts[0]++;
					break;
				case FIXEDVERTEX:
					fixedverts[0]++;
					continue;
				default:
					printf("bizarre vertex type %d in smooth, dying\n", kind);
					break;
				}

				/* try to smooth it */
				/* only record affected tets if we were given an output stack */
				worstqual = 1e-100;
				smoothed = smoothsinglevertex(currtpoint, threshold, worstqual, 
					smoothkinds, lapattempts, lapsuccesses, optattempts, optsuccesses);

				/* record number of successful smooths for this type of vertex */
				if (smoothed)
				{
					switch (kind)
					{
					case FREEVERTEX: 
						freeverts[1]++;
						break;
					case FACETVERTEX:
						facetverts[1]++;
						break;
					case SEGMENTVERTEX:
						segmentverts[1]++;
						break;
					case FIXEDVERTEX:
						fixedverts[1]++;
						break;
					}
				}
			}

			/* calculate the worst quality after smoothing */
			if (points != NULL)
				*minqualafter = qualitycal->minstackquality(tmesh, tetravec, qualitymetric);
			else
				*minqualafter = qualitycal->minmeshquality(tmesh, qualitymetric);

		} while ((*minqualafter - minqualbefore) > threshold);

		if (origqual < *minqualafter)
			return true;

		return false;
	}


	/* combination lap/opt smoothing for a single point ph */
	bool VertexSmoother::smoothsinglevertex(PointHandle ph, double threshold, double &worstout, int kinds,
		int &lapattempts, int &lapsuccesses, int &optattempts, int &optsuccesses)
	{
		bool noghosts = true;
		bool optsmoothed = false;
		bool lapsmoothed = false;
		double thisqual, startworstqual = HUGEFLOAT;
		Point tetp[4];

		/* a list of all the tets incident to this vertex */
		std::vector<TetraHandle> pointStar;
		/* a list of information about all the incident tets */
		std::vector<OptTet> incidenttets;
		OptTet opttet;

		/* check that the point exists */
		if (!tmesh->is_valid(ph))
			return false;

		tmesh->point_star(ph, pointStar);
		assert(pointStar.size() > 0);

		/* perform a pass of laplacian smoothing */
		++ lapattempts;
		lapsmoothed = lapsmooth(ph, threshold, worstout, kinds);
		/* if laplacian smoothing succeeds, return true */
		if (lapsmoothed)
		{
			++ lapsuccesses;
			return true;
		}

		/* if that fails to acheive desired quality, use non-smooth optimization-based smoothing */

		/* copy tags to incident tet data structure */
		for (unsigned int i = 0; i < pointStar.size(); i++)
		{
			/* copy vertex tags */
			opttet.handle = pointStar[i];
			incidenttets.push_back(opttet);

			/* compute starting worst quality */
			if (SMOOTHPARANOID)
			{
				tetPoints(tmesh, pointStar[i], tetp);
				thisqual = qualitycal->tetquality(tetp[0], tetp[1], tetp[2], tetp[3], qualitymetric);

				/* is this the worst quality incident tet? */
				if (thisqual < startworstqual) startworstqual = thisqual;
			}
		}

		//if (do_nonsmooth)
		//{
		//	++ optattempts;
		//	optsmoothed = nonsmooth(ph, incidenttets, worstout, kinds);

		//	if (optsmoothed) ++ optsuccesses;
		//}

		if (optsmoothed)
			return true;

		return false;
	}


	/* perform a pass of Laplacian smoothing. */
	bool VertexSmoother::lapsmooth(PointHandle ph, double threshold, double &worstout, int kinds)
	{
		Point tetpnt[4];
		Point newp;
		Point origpoint;       /* original point */
		SmoothVertex vinfo;    /* point information */
		SmoothVertex vinfotmp; /* temple point information */
		double origworst;      /* original worst quality */
		double afterworst;     /* after worst quality */
		double tetqual;
		std::vector<TetraHandle> pointStar;  /* tetras incident to the point */
		std::set<PointHandle> validpoint;  /* points incident to ph */
		std::set<PointHandle>::iterator validpoint_it;
		TetraMesh::HedronVertexIter hv_it;
		TetraMesh::HedronFaceIter hf_it;
		TetraMesh::HalfFaceVertexIter fv_it;
		TetraMesh::FaceHalfedgeIter fe_it;
		TetraMesh::HalfEdgeHandle heh;

		/* get point information */
		vertexinformation(vinfo, ph);

		/* check whether to try to smooth facet/segment/fixed vertices */
		if (vinfo.kind == FACETVERTEX && ((kinds & SMOOTHFACETVERTICES) == 0))
			return false;

		if (vinfo.kind == SEGMENTVERTEX && ((kinds & SMOOTHSEGMENTVERTICES) == 0))
			return false;

		/* save original point coordinate*/
		origpoint = tmesh->point(ph);

		/* get tetras incident to the point */
		tmesh->point_star(ph, pointStar);

		/* get valid incident points */
		validpoint.clear();

		switch (vinfo.kind)
		{
		case FREEVERTEX:
			/* all the points incident to the point are valid */
			for (unsigned int i = 0; i < pointStar.size(); ++i)
			{
				for (hv_it = tmesh->hedron_vertex_iter(pointStar[i]); hv_it; ++hv_it)
				{
					if (tmesh->point_handle(hv_it.handle()) != ph)
						validpoint.insert(tmesh->point_handle(hv_it.handle()));
				}
			}
			if (validpoint.size() < 3)
				return false;
			break;
		case FACETVERTEX:
			/* only the boundary points incident to the point are valid */
			for (unsigned int i = 0; i < pointStar.size(); ++i)
			{
				/* search each face in the tetra */
				for (hf_it = tmesh->hedron_face_iter(pointStar[i]); hf_it; ++hf_it)
				{
					/* if it is boundary, insert all its points */
					if (tmesh->is_boundary(hf_it.handle()) && tmesh->has_point(hf_it.handle(),ph))
					{
						fv_it = tmesh->half_face_vertex_iter(hf_it.handle());
						validpoint.insert(tmesh->point_handle(fv_it.handle())); ++fv_it;
						validpoint.insert(tmesh->point_handle(fv_it.handle())); ++fv_it;
						validpoint.insert(tmesh->point_handle(fv_it.handle()));
					}
				}
			}
			/* delete the smoothing point itself */
			validpoint_it = std::find(validpoint.begin(), validpoint.end(), ph);
			if(validpoint_it == validpoint.end())
				return false;
			validpoint.erase(validpoint_it);
			if (validpoint.size() < 3)
				return false;
			break;
		case SEGMENTVERTEX:
			/* only the segment points incident to the point are valid */
			for (unsigned int i = 0; i < pointStar.size(); ++i)
			{
				/* search each face in the tetra */
				for (hf_it = tmesh->hedron_face_iter(pointStar[i]); hf_it; ++hf_it)
				{
					/* if it is boundary, and contains the smoothing point */
					if (tmesh->is_boundary(hf_it.handle()) && tmesh->has_point(hf_it.handle(),ph))
					{
						/* iterate all the HalfEdge belong to the face */
						for (fe_it = tmesh->face_half_edge_iter(hf_it.handle()); fe_it; ++fe_it)
						{
							/* if the halfedge does not contain the smoothing point, go to next one */
							if (!tmesh->has_point(fe_it.handle(), ph))
								continue;
							/* fetch the opposite boundary half edge */
							heh = tmesh->mate_half_edge_handle(fe_it.handle());
							while (tmesh->has_radial_half_edge(heh))
							{
								heh = tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(heh));
							}
							/* if it is a sharp edge and contains the smoothing point */
							if (is_sharpedge(fe_it.handle(), heh))
							{
								validpoint.insert(tmesh->point_handle(tmesh->from_vertex_handle(fe_it.handle())));
								validpoint.insert(tmesh->point_handle(tmesh->to_vertex_handle(fe_it.handle())));
							}
						}
					}
				}
			}
			/* delete the smoothing point itself */
			validpoint_it = std::find(validpoint.begin(), validpoint.end(), ph);
			if (validpoint_it == validpoint.end())
				return false;
			validpoint.erase(validpoint_it);

			if (validpoint.size() != 2)
				return false;
			break;
		default:
			return false;
			break;
		}

		/* calculate original worst quality */
		origworst = qualitycal->minstackquality(tmesh, pointStar, qualitymetric);
		
		/* smooth point */
		newp[0] = newp[1] = newp[2] = 0.0;
		for (validpoint_it = validpoint.begin(); validpoint_it != validpoint.end(); ++ validpoint_it)
		{
			newp += tmesh->point(*validpoint_it);
		}
		newp /= validpoint.size();

		/* set the point to the new coordinate */
		tmesh->set_point(ph, newp);

		/* calculate worst quality after smoothing */
		afterworst = qualitycal->minstackquality(tmesh, pointStar, qualitymetric);
		worstout = afterworst;

		/*If the quality of tetras is improved, return true. */
		if (afterworst - origworst >= threshold)
		{
			/* record this vertex insertion in the journal */
			pointvec.clear();thvec.clear();phvec.clear();vhvec.clear();
			pointvec.push_back(origpoint);
			phvec.push_back(ph);
			journals->pushJournal(SMOOTHVERTEX, pointvec, thvec, phvec, vhvec);

			return true;
		}

		/* reset the point to the original place */
		tmesh->set_point(ph, origpoint);
		return false;
	}
}