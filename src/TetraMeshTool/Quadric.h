#ifndef TETRA_MESH_TOOL_QUADRIC_
#define TETRA_MESH_TOOL_QUADRIC_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/top.h>
#include <map>
namespace VolumeMesh
{

	/* surface error quadric */
	struct Quadric
	{
		double a2, ab, ac, ad;
		double     b2, bc, bd;
		double         c2, cd;
		double             d2;
		Point origpos;     /* original vertex position */
		int numfaces;            /* number of incident faces */
		double facesum;        /* sum of incident face areas */
		double edge2harm;      /* harmonic mean of squares of inc. edge lengths */
		bool hasquadric;
	};

	class QuadricContainer
	{
	public:
		QuadricContainer(): tmesh(NULL), qualitymetric(0),
			                 quadricoffset(0.8), quadricscale(300.0){}
		QuadricContainer(TetraMesh *tmesh_, int qualitymetric_): tmesh(tmesh_), qualitymetric(qualitymetric_),
			                                                     quadricoffset(0.8), quadricscale(300.0){}
		~QuadricContainer(){}
	public:
		void setMesh(TetraMesh *tmesh_)
		{
			tmesh = tmesh_;
		}

		void setQualityMetric(int qualitymetric_)
		{
			qualitymetric = qualitymetric_;
		}

		/* add a quadric for a specific vertex */
		void addquadric(PointHandle ph);

		/* get point quadric */
		Quadric quadric(PointHandle ph)
		{
			/* if we don't have the point's quadric, then add it firstly */
			if (!pointquadricmap.count(ph))
				addquadric(ph);
			return pointquadricmap[ph];
		}
		
		/* compute the quadric error at a vertex, normalized to
		be comparable to tetrahedron quality measures */
		double quadricerrortet(PointHandle ph);

		/* compute the quadric error for a query position of a vertex */
		double quadricerrorquery(PointHandle ph);

		/* compute the gradient of the quadric error, scaled for tet comparison */
		TetraMesh::Normal quadricgradtet(PointHandle ph);

		/* compute the gradient of the quadric error for query point */
		TetraMesh::Normal quadricgradquery(PointHandle ph);


	protected:
	private:
		/* Quadric factors */
		double quadricoffset;      /* quality to start every quadric at */
		double quadricscale;       /* factor to scale quadric by */

		TetraMesh *tmesh;
		std::map<PointHandle, Quadric> pointquadricmap;  /* point quadric map */
		int qualitymetric;
	};
}
#endif