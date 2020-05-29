#include <TetraMeshTool/NonsmoothQuality.h>
#include <TetraMeshTool/NonsmoothMath.h>
#include <math.h>
#include <fstream>
namespace VolumeMesh
{
	double NonsmoothQuality::getMinsine(NonsmoothPoint p[4])
	{
		NonsmoothPoint t, u, v;
		NonsmoothPoint facenormals[4];
		double anglecos[6];
		double anglesin[6];
		double minsin;
		t = p[1] - p[0];
		u = p[2] - p[0];
		v = p[3] - p[0];

		// get four normalized face normals
		facenormals[0] = unitVector((u - t) % (v - t));
		facenormals[1] = unitVector(v % u);
		facenormals[2] = unitVector(t % v);
		facenormals[3] = unitVector(u % t);

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
			if (minsin > anglesin[i])
				minsin = anglesin[i];
		}

		return minsin;
	}
	double NonsmoothQuality::getMinsine2(NonsmoothPoint p[4])
	{
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		NonsmoothPoint facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double sine2, minsine2;  /* the sine (squared) of the dihedral angle */
		int i, j, k, l;          /* loop indices */
		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = 6 * TetrahedronVolume(p);

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
				if ((facearea2[k] > 0) && (facearea2[l] > 0))
					sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
				else
				{
					sine2 = 0.0;
				}

				/* update minimum sine */
				if (sine2 < minsine2)
					minsine2 = sine2;
			}
		}
		double result = sqrt(minsine2) * pyrvolume;
		return result;
	}
	double NonsmoothQuality::getBiasedMinsine(NonsmoothPoint p[4])
	{
		NonsmoothPoint t, u, v;
		NonsmoothPoint facenormals[4];
		double anglecos[6];
		double anglesin[6];
		double minsin;
		t = p[1] - p[0];
		u = p[2] - p[0];
		v = p[3] - p[0];

		// get four normalized face normals
		facenormals[0] = unitVector((u - t) % (v - t));
		facenormals[1] = unitVector(v % u);
		facenormals[2] = unitVector(t % v);
		facenormals[3] = unitVector(u % t);

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
			// if biasis
			if (anglecos[i] < 0)
				anglesin[i] *= 0.7;
			if (minsin > anglesin[i])
				minsin = anglesin[i];
		}

		return minsin;

	}
	double NonsmoothQuality::getVomlength(NonsmoothPoint p[4])
	{
		double edgelength[3][4]; /* the lengths of each of the edges of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		int i, j, k, l;          /* loop indices */
		double edgelengthsum = 0.0;
		double lrms;             /* root mean squared of edge length */
		double qual;
		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = 6 * TetrahedronVolume(p);

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
		qual = ( pyrvolume) / (sqrt(2.0) * lrms * lrms * lrms);

		return qual;
	}
	double NonsmoothQuality::getZ(NonsmoothPoint p[4])
	{
		NonsmoothPoint ot,dt,ft;
		NonsmoothPoint crossdf, crossfo, crossod;
		NonsmoothPoint tettapex,tetorg,tetdest,tetfapex;
		NonsmoothPoint ct;
		double otlength, dtlength, ftlength;

		/* Use coordinates relative to the apex of the tetrahedron. */
		tettapex = p[0];
		tetorg = p[1];
		tetdest = p[2];
		tetfapex = p[3];
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
		ct.point_coord[0] = otlength * crossdf.point_coord[0] + dtlength * crossfo.point_coord[0] + ftlength * crossod.point_coord[0];
		ct.point_coord[1] = otlength * crossdf.point_coord[1] + dtlength * crossfo.point_coord[1] + ftlength * crossod.point_coord[1];
		ct.point_coord[2] = otlength * crossdf.point_coord[2] + dtlength * crossfo.point_coord[2] + ftlength * crossod.point_coord[2];
		/* Calculate the length of this vector, which is Z */
		return sqrt(ct | ct);
	}
	double NonsmoothQuality::getRadiusratio(NonsmoothPoint p[4])
	{
		NonsmoothPoint facenormal[4]; /* the normals of each face of the tet */
		double facearea2[4];     /* areas of the faces of the tet */
		double pyrvolume;        /* volume of tetrahedron */
		double Z;                /* quantity needed for circumradius */
		double facesum=0.0;       /* sum of the areas of the faces */
		double sign;
		double qual;

		/* calculate the volume*6 of the tetrahedron */
		pyrvolume = 6 * TetrahedronVolume(p);

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
		Z = getZ(p);

		/* now we are ready to compute the radius ratio, which is
		(108 * V^2) / Z (A1 + A2 + A3 + A4)
		(use 3 instead of 108 because pyrvolume = 6V)
		*/
		/* use sqrt for now... */
		sign = (pyrvolume < 0.0) ? -1.0 : 1.0;

		qual = sign * sqrt((3.0 * pyrvolume * pyrvolume) / (Z * facesum));
		return qual;
	}
}