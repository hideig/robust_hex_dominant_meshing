#include <TetraMeshTool/Quadric.h>
namespace VolumeMesh
{
	/* add a quadric for a newly inserted vertex */
	void QuadricContainer::addquadric(PointHandle ph)
	{
		Point newv, vptr[3];
		TetraMesh::Normal e1, e2, normal;
		double facearea, d, normfactor;
		std::vector<TetraMesh::HalfFaceHandle> faces;
		std::vector<TetraHandle> pointStar; /* incident tetras */
		TetraMesh::HedronFaceIter hf_it;
		TetraMesh::HalfFaceVertexIter fv_it;

		/* create quadric */
		Quadric *q = new Quadric;

		/* if the point is not boundary, it has no quadric*/
		if (!tmesh->is_boundary(ph))
		{
			q->hasquadric = false;
		}
		else
		{
			/* initialize quadric */
			q->hasquadric = true;

			q->a2 = q->ab = q->ac = q->ad 
				= q->b2 = q->bc = q->bd 
				= q->c2 = q->cd 
				= q->d2 = 0.0;

			q->numfaces = 0;
			q->facesum = 0.0;
			q->edge2harm = 0.0;

			newv = tmesh->point(ph);
			q->origpos = newv;

			/* fetch incident boundary faces */
			tmesh->point_star(ph, pointStar);
			for (unsigned int i = 0; i < pointStar.size(); ++i)
			{
				for (hf_it = tmesh->hedron_face_iter(pointStar[i]); hf_it; ++hf_it)
				{
					if (tmesh->is_boundary(hf_it.handle()) && tmesh->has_point(hf_it.handle(), ph))
					{
						faces.push_back(hf_it.handle());
					}
				}
			}

			/* accumulate fundamental quadric for each incident face */
			for (unsigned int i = 0; i < faces.size(); ++i)
			{
				/* get the actual vertices */
				fv_it = tmesh->half_face_vertex_iter(faces[i]);
				vptr[0] = tmesh->point(fv_it.handle()); ++fv_it;
				vptr[1] = tmesh->point(fv_it.handle()); ++fv_it;
				vptr[2] = tmesh->point(fv_it.handle());

				/* compute face normal */
				e1 = (vptr[1] - vptr[0]);
				e2 = (vptr[2] - vptr[0]);
				normal =  e1 % e2;

				/* face area is 1/2 the length of the cross product */
				facearea = (normal | normal) / 2.0;

				/* normalize the normal */
				normal.normalize();

				/* compute the orthogonal distance from the plane of this face to
				the origin */
				d = -(normal | vptr[0]);

				q->numfaces++;
				q->facesum += facearea;
				/* add on 1/2 of 1 / l^2, because every edge will be counted twice */
				if (e1.norm() > 0.0 && e2.norm() > 0.0)
				{
					q->edge2harm += (0.5 / (e1|e1));
					q->edge2harm += (0.5 / (e2|e2));
				}

				/* accumulate the fundamental quadric from this face */
				/* normal = [a b c] */
				q->a2 += normal[0] * normal[0] * facearea;
				q->ab += normal[0] * normal[1] * facearea;
				q->ac += normal[0] * normal[2] * facearea;
				q->ad += normal[0] * d * facearea;
				q->b2 += normal[1] * normal[1] * facearea;
				q->bc += normal[1] * normal[2] * facearea;
				q->bd += normal[1] * d * facearea;
				q->c2 += normal[2] * normal[2] * facearea;
				q->cd += normal[2] * d * facearea;
				q->d2 += d * d * facearea;
			}

			/* compute normalization */
			/* quadric must have at least 3 surrounding faces */
			assert(q->numfaces >= 3);

			/* compute harmonic mean */
			q->edge2harm = ((double) q->numfaces) / q->edge2harm;

			/* compute normalization factor */
			normfactor = q->edge2harm * q->facesum;

			/* if the facesum is zero, bail */
			if (q->facesum <= 0.0)
			{
				q->hasquadric = false;
			}

			else
			{
				assert(q->edge2harm != 0.0);
				assert(q->facesum != 0.0);
				assert(normfactor != 0.0);

				/* scale quadric by normalization factor */
				q->a2 /= normfactor;
				q->ab /= normfactor;
				q->ac /= normfactor;
				q->ad /= normfactor;
				q->b2 /= normfactor;
				q->bc /= normfactor;
				q->bd /= normfactor;
				q->c2 /= normfactor;
				q->cd /= normfactor;
				q->d2 /= normfactor;
			}
		}
		pointquadricmap[ph] = *q;
	}


	/* compute the quadric error at a vertex, normalized to
	be comparable to tetrahedron quality measures */
	double QuadricContainer::quadricerrortet(PointHandle ph)
	{
		if (!pointquadricmap.count(ph))
			addquadric(ph);

		/* return the quadric error of this vertex scaled
		by the quadric scale and offset by the quadric offset */
		double qe = quadricoffset - (quadricscale * quadricerrorquery(ph));

		/* scale down so that perfect corresponds to equilateral */
		if (qualitymetric == QUAL_MINSINE || qualitymetric == QUAL_WARPEDMINSINE)
		{
			qe *= SINEEQUILATERAL;
		}

		/* don't return negative qualities */
		if (qe < 0.0) return 0.0;
		return qe;
	}

	/* compute the quadric error for a query position of a vertex */
	double QuadricContainer::quadricerrorquery(PointHandle ph)
	{
		if (!pointquadricmap.count(ph))
			addquadric(ph);

		Quadric *q;
		double quad;
		Point v;

		/* fetch the quadric for this vertex */
		q = new Quadric;
		assert(q != NULL);

		/* return zero if this is not a surface vertex */
		q = &quadric(ph);
		if (q->hasquadric == false) return 0.0;

		/* evaluate quadric */
		/* Q(v) = v'Av + 2b'v + c */
		v = tmesh->point(ph);
		quad = v[0]*v[0]*q->a2 + 2*v[0]*v[1]*q->ab + 2*v[0]*v[2]*q->ac + 2*v[0]*q->ad
			 + v[1]*v[1]*q->b2 + 2*v[1]*v[2]*q->bc + 2*v[1]*q->bd
		     + v[2]*v[2]*q->c2 + 2*v[2]*q->cd
			 + q->d2;

		return quad;
	}

	/* compute the gradient of the quadric error, scaled for tet comparison */
	TetraMesh::Normal QuadricContainer::quadricgradtet(PointHandle ph)
	{
		if (!pointquadricmap.count(ph))
			addquadric(ph);
		return quadricgradquery(ph);
	}

	/* compute the gradient of the quadric error for query point */
	TetraMesh::Normal QuadricContainer::quadricgradquery(PointHandle ph)
	{
		Point v = tmesh->point(ph);
		Quadric *q;
		TetraMesh::Normal grad;

		/* fetch the quadric for this vertex */
		q = new Quadric;
		assert(q != NULL);

		/* return a zero gradient if this is not a surface vertex */
		if (q->hasquadric == false)
		{
			grad[0] = grad[1] = grad[2] = 0.0;
			return grad;
		}

		/* grad(Q) = 2Av + 2b */
		/* A = nn', b = dn */
		/* so */
		/* grad(Q) = 2 [xa2 + yab + zac + ad] */
		/*             [xab + yb2 + zbc + bd] */
		/*             [xac + ybc + zc2 + cd] */
		/* negative because we want to *reduce* quadric error... */
		grad[0] = -2.0 * (v[0]*q->a2 + v[1]*q->ab + v[2]*q->ac + q->ad);
		grad[1] = -2.0 * (v[0]*q->ab + v[1]*q->b2 + v[2]*q->bc + q->bd);
		grad[2] = -2.0 * (v[0]*q->ac + v[1]*q->bc + v[2]*q->c2 + q->cd);

		return grad;
	}

}