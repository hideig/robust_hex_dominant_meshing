#include <TetraMeshTool/QualityCalculator.h>
#include <algorithm>
#include <TetraMeshTool/StarBase.h>


namespace VolumeMesh
{	

	/* compute the (square) of the minimum sine
	 * of all the dihedral angles in the tet defined
	 * by the four vertices (p1, p2, p3, p4).
	 * The order of vertices should be consistent with VolumeMesh defined.
	*/
	//double QualityCalculator::minsine(Point p1, Point p2, Point p3, Point p4)
	//{
	//	TetraMesh::Normal t, u, v;
	//	TetraMesh::Normal facenormals[4];
	//	double anglecos[6];
	//	double anglesin[6];
	//	double minsin;

	//	t = p2 - p1;
	//	u = p3 - p1;
	//	v = p4 - p1;

	//	// get four normalized face normals
	//	facenormals[0] = ((u - t) % (v - t)).normalize();
	//	facenormals[1] = (v % u).normalize();
	//	facenormals[2] = (t % v).normalize();
	//	facenormals[3] = (u % t).normalize();

	//	// get six tetra dihedral angles cosine
	//	anglecos[0] = facenormals[0] | facenormals[1];
	//	anglecos[1] = facenormals[0] | facenormals[2];
	//	anglecos[2] = facenormals[0] | facenormals[3];
	//	anglecos[3] = facenormals[1] | facenormals[2];
	//	anglecos[4] = facenormals[1] | facenormals[3];
	//	anglecos[5] = facenormals[2] | facenormals[3];

	//	// get six tetra dihedral angles sine and the smallest one
	//	minsin = 1.0;
	//	for (int i = 0; i < 6; i ++)
	//	{
	//		anglesin[i] = pow((1-anglecos[i]*anglecos[i]), 0.5) * (3 / 2 / pow(2, 0.5));
	//		if (minsin > anglesin[i])
	//			minsin = anglesin[i];
	//	}

	//	return minsin;
	//}

	double QualityCalculator::minsine(Point p1, Point p2, Point p3, Point p4)
	{
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		TetraMesh::Normal facenormal[4]; /* the normals of each face of the tet */
		Point p[4];                      /* the four points of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double sine2, minsine2;  /* the sine (squared) of the dihedral angle */
		int i, j, k, l;          /* loop indices */
		
		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0]) % (p[2]-p[0])) | (p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
		{
			return 0.0;
		}

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			facearea2[i] = facenormal[i] | facenormal[i];
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
		}

		minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */

		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) 
		{
			for (j = i + 1; j < 4; j++) 
			{
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the expression for minimum sine, squared, over 4 
				The reason it's over 4 is because the area values we have
				are actually twice the area squared */
				/* if either face area is zero, the sine is zero */
				if (facearea2[k] > 0 && facearea2[l] > 0)
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				else
					sine2 = 0.0;

				/* update minimum sine */
				if (sine2 < minsine2)
					minsine2 = sine2;
			}
		}
		return sqrt(minsine2) * pyrvolume;
	}

	double QualityCalculator::biasedminsine(Point p1, Point p2, Point p3, Point p4)
	{
		TetraMesh::Normal t, u, v;
		TetraMesh::Normal facenormals[4];
		double anglecos[6];
		double anglesin[6];
		double minsin;

		t = p2 - p1;
		u = p3 - p1;
		v = p4 - p1;

		// get four normalized face normals
		facenormals[0] = ((u - t) % (v - t)).normalize();
		facenormals[1] = (v % u).normalize();
		facenormals[2] = (t % v).normalize();
		facenormals[3] = (u % t).normalize();

		// get six tetra dihedral angles cosine
		anglecos[0] = facenormals[0] | facenormals[1];
		anglecos[1] = facenormals[0] | facenormals[2];
		anglecos[2] = facenormals[0] | facenormals[3];
		anglecos[3] = facenormals[1] | facenormals[2];
		anglecos[4] = facenormals[1] | facenormals[3];
		anglecos[5] = facenormals[2] | facenormals[3];

		// get six tetra dihedral angles sine and the smallest one
		minsin = 1.0;
		for (int i = 0; i < 6; i ++)
		{
			anglesin[i] = pow((1-anglecos[i]*anglecos[i]), 0.5) * (3 / 2 / pow(2, 0.5));
			// if biase
			if (anglecos[i] < 0)
				anglesin[i] *= 0.7;
			if (minsin > anglesin[i])
				minsin = anglesin[i];
		}

		return minsin;
	}

	/* compute the mean of the sines
	 * of all the dihedral angles in the tet defined
	 * by the four vertices (p1, p2, p3, p4).
	 * The order of vertices should be consistent with VolumeMesh defined.
	*/
	double QualityCalculator::meansine(Point p1, Point p2, Point p3, Point p4)
	{
		Point p[4];      /* the vertices of the tet */
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		Vec3d facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double sine2;            /* the sine (squared) of the dihedral angle */
		double sinesum=0.0;      /* the accumulating sum of the sines */
		int i, j, k, l;          /* loop indices */

		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
			return 0.0;

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			facearea2[i] = facenormal[i] | facenormal[i];
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
		}

		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the expression for minimum sine, squared, over 4 
				The reason it's over 4 is because the area values we have
				are actually twice the area squared */
				/* if either face area is zero, the sine is zero */
				if (facearea2[k] > 0 && facearea2[l] > 0)
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				else
					sine2 = 0.0;

				/* accumulate sine */
				sinesum += sqrt(sine2);
			}
		}

		/* average sine */
		return (sinesum / 6.0) * pyrvolume;
	}

	/* compute the (square) of the minimum sine
	of all the dihedral angles in the tet defined
	by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	double QualityCalculator::minsineandedgeratio(Point p1, Point p2, Point p3, Point p4)
	{
		Point p[4];      /* the vertices of the tet */
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		Vec3d facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double sine2,minsine2;   /* the sine (squared) of the dihedral angle */
		double sinesum=0.0;      /* the accumulating sum of the sines */
		double minsine;
		double shortest = HUGEFLOAT;
		double longest = 0.0;     /* shortest and longest edges in teh tet */
		double edgeratio;         /* ratio of shortest to longest edge */
		int i, j, k, l;          /* loop indices */

		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
			return 0.0;

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			facearea2[i] = facenormal[i] | facenormal[i];
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
			{
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
				/* keep track of longest and shortest edge */
				if (edgelength[i][j] > longest) longest = edgelength[i][j];
				if (edgelength[i][j] < shortest) shortest = edgelength[i][j];
			}
		}

		minsine2 =  10.0e10;     /* start with absurdly big value for sine */


		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the expression for minimum sine, squared, over 4 
				The reason it's over 4 is because the area values we have
				are actually twice the area squared */
				/* if either face area is zero, the sine is zero */
				if (facearea2[k] > 0 && facearea2[l] > 0)
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				else
					sine2 = 0.0;
				/* update minimum sine */
				if (sine2 < minsine2)
				{
					minsine2 = sine2;
				}
			}
		}

		/* edge ratio, scaled down for parity with sin of equilateral tet's dihedrals */
		edgeratio = sqrt(shortest / longest) * SINEEQUILATERAL;
		minsine = sqrt(minsine2) * pyrvolume;

		if (edgeratio < minsine) 
			return edgeratio;

		return minsine;
	}

	/* compute Z, a quantity associated with circumradius computation
	TODO this code is lifted from Jonathan's tetcircumcenter computation
	in primitives.c */
	double QualityCalculator::getZ(Point tetorg, Point tetdest, Point tetfapex, Point tettapex)
	{
		TetraMesh::Normal ot,dt,ft;
		TetraMesh::Normal crossdf, crossfo, crossod;
		TetraMesh::Normal ct;
		double otlength, dtlength, ftlength;

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
		ct[0] = otlength * crossdf[0] + dtlength * crossfo[0] + ftlength * crossod[0];
		ct[1] = otlength * crossdf[1] + dtlength * crossfo[1] + ftlength * crossod[1];
		ct[2] = otlength * crossdf[2] + dtlength * crossfo[2] + ftlength * crossod[2];
		/* Calculate the length of this vector, which is Z */
		return sqrt(ct | ct);
	}

	/* the inradius to circumradius ratio */
	double QualityCalculator::radiusratio(Point p1, Point p2, Point p3, Point p4)
	{
		Point p[4];      /* the vertices of the tet */
		Vec3d facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double Z;                /* quantity needed for circumradius */
		double facesum=0.0;       /* sum of the areas of the faces */
		double sign;
		double qual;

		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
			return 0.0;

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (int i = 0; i < 4; ++i)
		{
			facearea2[i] = facenormal[i] | facenormal[i];
			facesum += sqrt(facearea2[i]) * 0.5;
		}

		/* compute Z */
		Z = getZ(p[0], p[1], p[2], p[3]);

		/* now we are ready to compute the radius ratio, which is
		(108 * V^2) / Z (A1 + A2 + A3 + A4)
		(use 3 instead of 108 because pyrvolume = 6V)
		*/
		/* use sqrt for now... */
		sign = (pyrvolume < 0.0) ? -1.0 : 1.0;

		qual = sign * sqrt((3.0 * pyrvolume * pyrvolume) / (Z * facesum));
		return qual;
	}

	/* compute the ratio of the tet volume to the cube of
	the rms edge length */
	double QualityCalculator::vlrms3ratio(Point p1, Point p2, Point p3, Point p4)
	{
		Point p[4];      /* the vertices of the tet */
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		int i, j, k, l;          /* loop indices */
		double edgelengthsum = 0.0;
		double lrms;             /* root mean squared of edge length */
		double qual;

		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
			return 0.0;

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
		}


		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				edgelengthsum += edgelength[i][j];
			}
		}

		/* compute the root mean square */
		lrms = sqrt((1.0 / 6.0) * edgelengthsum);

		/* compute the normalized ratio of volume to lrms^3 */
		qual = (sqrt(2.0) * pyrvolume) / (lrms * lrms * lrms);

		return qual;
	}

	/* compute the (square) of the minimum sine
	of all the dihedral angles in the tet defined
	by the four vertices (vtx1, vtx2, vtx3, vtx4)
	*/
	double QualityCalculator::warpedminsine(Point p1, Point p2, Point p3, Point p4)
	{
		Point p[4];      /* the vertices of the tet */
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		TetraMesh::Normal facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double sine2, minsine2;  /* the sine (squared) of the dihedral angle */
		int i, j, k, l;          /* loop indices */

		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
		{
			return 0.0;
		}

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			facearea2[i] = facenormal[i] | facenormal[i];
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
		}

		minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */


		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the expression for minimum sine, squared, over 4 
				The reason it's over 4 is because the area values we have
				are actually twice the area squared */
				/* if either face area is zero, the sine is zero */
				if (facearea2[k] > 0 && facearea2[l] > 0)
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				else
					sine2 = 0.0;

				/* check whether this angle is obtuse */
				if ((facenormal[k] | facenormal[l]) > 0)
				{
					/* if so, warp it down */
					sine2 = (sinewarpfactor) * sqrt(sine2);
					sine2 *= sine2;
				}

				/* update minimum sine */
				if (sine2 < minsine2)
				{
					minsine2 = sine2;
				}
			}
		}

		return sqrt(minsine2) * pyrvolume;
	}

	/* compute the minimum or maximum angle of the tet defined
	 * by the four vertices (vtx1, vtx2, vtx3, vtx4)
	 */
	double QualityCalculator::minmaxangle(Point p1, Point p2, Point p3, Point p4, bool max)
	{
		Point p[4];      /* tet vertices */
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		TetraMesh::Normal facenormal[4]; /* the normals of each face of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		int i, j, k, l;          /* loop indices */
		double minangle = HUGEFLOAT;
		double maxangle = 0.0;
		double angle, tantheta;
		double dotproduct;


		p[0] = p1;p[1] = p2;p[2] = p3;p[3] = p4;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = ((p[1]-p[0])%(p[2]-p[0]))|(p[3]-p[0]);

		/* if the volume is zero, the quality is zero, no reason to continue */
		if (pyrvolume == 0.0)
		{
			return 0.0;
		}

		/* compute the normal for each face, face(i) is opposite to vertex(i) */
		facenormal[0] = (p[3]-p[1]) % (p[2]-p[1]);
		facenormal[1] = (p[3]-p[0]) % (p[2]-p[0]);
		facenormal[2] = (p[1]-p[0]) % (p[3]-p[0]);
		facenormal[3] = (p[2]-p[0]) % (p[1]-p[0]);

		/* compute (2 *area)^2 for this face */
		for (i = 0; i < 4; ++i)
		{
			/* compute edge lengths (squared) */
			for (j = i + 1; j < 4; j++)
				edgelength[i][j] = (p[i]-p[j]) | (p[i]-p[j]);
		}


		/* for each edge in the tetrahedron */
		for (i = 0; i < 3; i++) {
			for (j = i + 1; j < 4; j++) {
				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
				l = 6 - i - j - k;

				/* compute the tangent of the angle using the tangent formula:

				tan(theta_ij) = - 6 * V * l_ij
				--------------
				dot(n_k, n_l)

				because this formula is accurate in the entire range.
				*/
				dotproduct = facenormal[k] | facenormal[l];

				if (dotproduct != 0.0)
				{
					tantheta = (-pyrvolume * sqrt(edgelength[i][j])) / dotproduct;

					/* now compute the actual angle */
					angle = atan(tantheta);
				}
				else
				{
					angle = PI / 2.0;
				}

				/* adjust angle for sign of dot product */
				if (dotproduct > 0)
				{
					angle += PI;
				}

				/* make negative angles positive */
				if (angle < 0)
				{
					angle += 2.0 * PI;
				}

				if (dotproduct == 0.0) angle = PI / 2.0;

				if (angle < minangle) minangle = angle;
				if (angle > maxangle) maxangle = angle;
			}
		}

		/*
		assert(radtodeg(maxangle) <= 180.0);
		assert(minangle >= 0.0);
		*/
		if (max) return radtodeg(maxangle);
		return radtodeg(minangle);

	}

	double QualityCalculator::tetquality(Point p1, Point p2, Point p3, Point p4, int qualityMetric)
	{
		double quality = 0.0; /* the quality of this tetrahedron */
		switch (qualityMetric)
		{
		case QUAL_MINSINE:
			quality = minsine(p1, p2, p3, p4);
			break;
		case QUAL_BIASEDMINSINE:
			quality = biasedminsine(p1, p2, p3, p4);
			break;
		case QUAL_MEANSINE:
			quality = meansine(p1, p2, p3, p4);
			break;
		case  QUAL_MINSINEANDEDGERATIO:
			quality = minsineandedgeratio(p1, p2, p3, p4);
			break;
		case  QUAL_RADIUSRATIO:
			quality = radiusratio(p1, p2, p3, p4);
			break;
		case QUAL_VLRMS3RATIO:
			quality = vlrms3ratio(p1, p2, p3, p4);
			break;
		case QUAL_WARPEDMINSINE:
			quality = warpedminsine(p1, p2, p3, p4);
			break;
		case QUAL_MINANGLE:
			quality = minmaxangle(p1, p2, p3, p4, false);
			break;
		case QUAL_MAXANGLE:
			quality = minmaxangle(p1, p2, p3, p4, true);
			break;
		default:
			printf("Undefined quality measure %d !\n", qualityMetric);
			double(1);
		}
		return quality;
	}

	/* compute the mean and minimum element
	 * qualities in the mesh (multiple thresholded means)
	 */
	void QualityCalculator::meshquality(TetraMesh *mesh, int qualmeasure, double *minqual, TetraStack *tetstack)
	{
		/* fill the stack of tets with all tets in the mesh */
		fillstackqual(mesh, qualmeasure, minqual, HUGEFLOAT, tetstack);
	}

	/* given a mesh and a quality threshold,
	 * return a stack of tets with quality 
	 * at or below that threshold. this function
	 * assumes that the stack has already been 
	 * initialized.
	 */
	void QualityCalculator::fillstackqual(TetraMesh *mesh, int qualmeasure, double *minqual, 
		                                  double threshold, TetraStack *tetstack)
	{
		TetraMesh::HedronIter h_iter;  /*the iterator of hedrons in mesh*/
		double quality;                /* the quality of the current tet */
		ImproveTetra itet;             /*the stack tetra*/
		int numtets = 0;
		Point tetpoint[4];

		/* make sure the stack is empty */
		if (tetstack != NULL)
			tetstack->clear();

		*minqual = HUGEFLOAT;

		 /* for each tetrahedron in mesh */
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++h_iter)
		{
			if (!tmesh->is_valid(h_iter.handle()))
				continue;
			/* get the points of current tetrahedron*/
			TetraMesh::HedronVertexIter hv_iter = mesh->hedron_vertex_iter(h_iter);
			tetpoint[0] = mesh->point(hv_iter.handle()); ++hv_iter;
			tetpoint[1] = mesh->point(hv_iter.handle()); ++hv_iter;
			tetpoint[2] = mesh->point(hv_iter.handle()); ++hv_iter;
			tetpoint[3] = mesh->point(hv_iter.handle()); 

 			/* compute the quality of this tet */
			quality = tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualmeasure);

			/* keep track of minimum quality */
			if (quality < *minqual) *minqual = quality;

			/* if this tet is worse than the threshold */
			if (tetstack != NULL && quality < threshold)
			{
				/* track thresholded mean qualities only for tets actually included in the stack */
				++numtets;

				/* push this tet on the "to be fixed" stack */
				itet.handle = h_iter.handle();
				itet.quality = quality;
				tetstack->push_back(itet);
			}
		}
	}

	/* find the meand and minimum qualities
	 * of all the tets in the stack (that still exist) */
	void QualityCalculator::stackquality(TetraMesh *mesh, TetraStack &tetstack, int qualmeasure, double meanqual[], double *minqual)
	{
		ImproveTetra *tet;   /* current tet */
		int nonexist = 0;
		double worstqual = HUGEFLOAT;
		int numtets = 0;
		double newqual;

		for (int j = 0; j < NUMMEANTHRESHOLDS; ++j)
			meanqual[j] = 0.0;

		for (unsigned int i = 0; i < tetstack.size(); ++i)
		{
			/* fetch this item from the stack */
			tet = &tetstack[i];

			/* check that it exists */
			if (!mesh->is_valid(tet->handle))
			{
				nonexist++;
				continue;
			}

			/* it exists, compute it's quality */
			TetraMesh::HedronVertexIter hv_iter = mesh->hedron_vertex_iter(tet->handle);
			Point tetpoint[4];
			tetpoint[0] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[1] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[2] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[3] = mesh->point(hv_iter.handle());
 			newqual = tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualmeasure);

			/* track thresholded mean qualities only for tets actually included in the stack */
			numtets++;
			for (int j=0; j<NUMMEANTHRESHOLDS; j++)
			{
				meanqual[j] += (newqual < meanthresholds[qualmetric][j]) ? newqual : meanthresholds[qualmetric][j];
			}

			/* is this a new low ? */
			if (newqual < worstqual) worstqual = newqual;

		}

		*minqual = worstqual;

		/* compute thresholded means */
		for (int i=0; i<NUMMEANTHRESHOLDS; i++)
		{
			meanqual[i] /= (double) numtets;
		}
	}

	void QualityCalculator::stackquality(TetraMesh *mesh, std::vector<TetraHandle> &tetstack, 
		                                 int qualmeasure, double meanqual[], double *minqual)
	{
		TetraHandle tet;   /* current tet */
		int nonexist = 0;
		double worstqual = HUGEFLOAT;
		int numtets = 0;
		double newqual;

		for (int j = 0; j < NUMMEANTHRESHOLDS; ++j)
			meanqual[j] = 0.0;

		for (unsigned int i = 0; i < tetstack.size(); ++i)
		{
			/* fetch this item from the stack */
			tet = tetstack[i];

			/* check that it exists */
			if (!mesh->is_valid(tet))
			{
				nonexist++;
				continue;
			}

			/* it exists, compute it's quality */
			TetraMesh::HedronVertexIter hv_iter = mesh->hedron_vertex_iter(tet);
			Point tetpoint[4];
			tetpoint[0] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[1] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[2] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[3] = mesh->point(hv_iter.handle());
			newqual = tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualmeasure);

			/* track thresholded mean qualities only for tets actually included in the stack */
			numtets++;
			for (int j=0; j<NUMMEANTHRESHOLDS; j++)
			{
				meanqual[j] += (newqual < meanthresholds[qualmetric][j]) ? newqual : meanthresholds[qualmetric][j];
			}

			/* is this a new low ? */
			if (newqual < worstqual) worstqual = newqual;

		}

		*minqual = worstqual;

		/* compute thresholded means */
		for (int i=0; i<NUMMEANTHRESHOLDS; i++)
		{
			meanqual[i] /= (double) numtets;
		}
	}

	double QualityCalculator::minstackquality(TetraMesh *mesh, const std::vector<TetraHandle> &tetstack, int qualmeasure)
	{
		TetraHandle tet;   /* current tet */
		double worstqual = HUGEFLOAT;
		double newqual;

		for (unsigned int i = 0; i < tetstack.size(); ++i)
		{
			/* fetch this item from the stack */
			tet = tetstack[i];

			/* check that it exists */
			if (!mesh->is_valid(tet))
				continue;

			/* it exists, compute it's quality */
			TetraMesh::HedronVertexIter hv_iter = mesh->hedron_vertex_iter(tet);
			Point tetpoint[4];
			tetpoint[0] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[1] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[2] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[3] = mesh->point(hv_iter.handle());
			newqual = tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualmeasure);

			/* is this a new low ? */
			if (newqual < worstqual) worstqual = newqual;
		}

		return worstqual;
	}


	double QualityCalculator::minstackquality(TetraMesh *mesh, TetraStack &tetstack, int qualmeasure)
	{
		TetraHandle tet;   /* current tet */
		double worstqual = HUGEFLOAT;
		double newqual;

		for (unsigned int i = 0; i < tetstack.size(); ++i)
		{
			/* fetch this item from the stack */
			tet = tetstack[i].handle;

			/* check that it exists */
			if (!mesh->is_valid(tet))
				continue;

			/* it exists, compute it's quality */
			TetraMesh::HedronVertexIter hv_iter = mesh->hedron_vertex_iter(tet);
			Point tetpoint[4];
			tetpoint[0] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[1] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[2] = mesh->point(hv_iter.handle());++ hv_iter;
			tetpoint[3] = mesh->point(hv_iter.handle());
			newqual = tetquality(tetpoint[0], tetpoint[1], tetpoint[2], tetpoint[3], qualmeasure);
			tetstack[i].quality = newqual;

			/* is this a new low ? */
			if (newqual < worstqual) worstqual = newqual;
		}

		return worstqual;
	}

	/* given an input mesh, calculate the worst quality value */
	double QualityCalculator::minmeshquality(TetraMesh *mesh, int qualmetric)
	{
		double qual, minqual;
		Point p[4];
		TetraMesh::HedronIter h_iter;
		minqual = 1.0;
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		{
			if (!tmesh->is_valid(h_iter.handle()))
				continue;
			tetPoints(mesh, h_iter.handle(), p);
			qual = tetquality(p[0], p[1], p[2], p[3], qualmetric);
			if (qual < minqual)
				minqual = qual;
		}
		return minqual;
	}

	/* given an input mesh, find the worst "input" angle.
	that is, find the smallest angle between two faces
	of the boundary */
	double QualityCalculator::worstinputangle(TetraMesh *mesh)
	{
		double angle;                     /* angle between two surface faces */
		double worstangle = 2.0 * PI;     /* worst input angle */
		//double minqual, meanqual[NUMMEANTHRESHOLDS];         /* quality of the worst tet in the mesh */
		int edgelist[6][2];               /* list of boundary edges for a tet */
		int edgefaces[6][2];
		TetraMesh::HedronIter h_iter;

		// set edge vertex handle
		edgelist[0][0] = 0;edgelist[0][1] = 1;
		edgelist[1][0] = 0;edgelist[1][1] = 2;
		edgelist[2][0] = 0;edgelist[2][1] = 3;
		edgelist[3][0] = 1;edgelist[3][1] = 2;
		edgelist[4][0] = 1;edgelist[4][1] = 3;
		edgelist[5][0] = 2;edgelist[5][1] = 3;

		// set face incident to the edge
		edgefaces[0][0] = 2;edgefaces[0][1] = 3;
		edgefaces[1][0] = 1;edgefaces[1][1] = 3;
		edgefaces[2][0] = 1;edgefaces[2][1] = 2;
		edgefaces[3][0] = 0;edgefaces[3][1] = 3;
		edgefaces[4][0] = 0;edgefaces[4][1] = 2;
		edgefaces[5][0] = 0;edgefaces[5][1] = 1;

		/* go through each tet on the stack */
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		{
			/* for each boundary edge */
			for (int i = 0; i < 6; i++)
			{
				/* compute the angle between the boundary faces */
				angle = getboundaryedgeangle(mesh, h_iter.handle(), edgelist[i][0], edgelist[i][1], edgefaces[i][0], edgefaces[i][1]);

				/* if this angle is smaller than what we've seen, update */
				if (angle < worstangle)
				{
					worstangle = angle;
				}
			}
		}
		return worstangle; 
	}
	
	/* given the two vertices of a known boundary edge,
	compute the angle between it's two incident
	boundary faces */
	double QualityCalculator::getboundaryedgeangle(TetraMesh *mesh, TetraHandle th, int v1, int v2, int vleft, int vright)
	{
		Point point[4];              /* the actual positions of the four vertices*/
		TetraMesh::Normal e[3];      /* edges for normal computation */
		TetraMesh::Normal norm[2];   /* face normals */
		double edgelength;       /* length of boundary edge */
		double volume;           /* volume of tet formed by four vertices */
		double dotprod;          /* the dot product of the two inward face normals */
		double tantheta;         /* tangent of the angle */
		double angle;            /* angle we are computing */

		/* fetch actual vertex values */
		point[0] = mesh->point(VertexHandle(th.idx()*4+v1));
		point[1] = mesh->point(VertexHandle(th.idx()*4+v2));
		point[2] = mesh->point(VertexHandle(th.idx()*4+vleft));
		point[3] = mesh->point(VertexHandle(th.idx()*4+vright));

		/* e[0] from v1 to v2 */
		e[0] = point[1] - point[0];
		/* e[1] from v1 to right */
		e[1] = point[2] - point[0];
		/* e[2] from v1 to left */
		e[2] = point[3] - point[0];

		/* compute the length of the edge in question */
		edgelength = e[0].norm();

		/* compute face normals, pointing out of the mesh interior */
		norm[0] = e[1] % e[0];
		norm[1] = e[0] % e[2];

		/* compute 6 * volume of the tetrahedron with these four vertices */
		volume = (e[0] % e[1]) | e[2];

		/* find the dot product of the two vectors */
		dotprod = norm[0] | norm[1];

		/* if the dot product is zero, then the angle is +-90 degrees */
		if (dotprod == 0)
		{
			/* if the volume is negative, angle is 270 degrees */
			angle = (volume > 0) ? PI / 2.0 : (3.0 * PI) / 2.0;
		}
		else
		{
			/* compute the tangent of the angle using the tangent formula:

			tan(theta_ij) = - 6 * V * l_ij
			--------------
			dot(n_k, n_l)

			because this formula is accurate in the entire range.
			*/
			tantheta = (-volume * edgelength) / dotprod;

			/* now compute the actual angle */
			angle = atan(tantheta);

			/* adjust angle for sign of dot product */
			if (dotprod > 0)
			{
				angle += PI;
			}

			/* make negative angles positive */
			if (angle < 0)
			{
				angle += 2.0 * PI;
			}

			/* zero angles must actually be 360 degrees...?? */
			if (angle == 0)
			{
				angle = 2.0 * PI;
			}
		}
		return angle;
	}

	/* return the worst quality of all elements in the mesh */
	double QualityCalculator::worstquality(TetraMesh *mesh)
	{
		double worstqual;
		double quality;
		TetraMesh::HedronIter h_iter;
		Point p[4];
		worstqual = 1.0;
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		{
			tetPoints(mesh, h_iter.handle(), p);
			quality = tetquality(p[0], p[1], p[2], p[3], qualmetric);
			if (quality < worstqual)
				worstqual = quality;
		}
		return worstqual;
	}

	/* return the worst quality of all elements in the vector */
	double QualityCalculator::worstquality(const std::vector<TetraHandle> &tetvec)
	{
		double worstqual;
		double quality;
		Point p[4];
		worstqual = 1.0;
		for (unsigned int i = 0; i < tetvec.size(); i++)
		{
			tetPoints(tmesh, tetvec[i], p);
			quality = tetquality(p[0], p[1], p[2], p[3], qualmetric);
			if (quality < worstqual)
				worstqual = quality;
		}
		return worstqual;
	}
}