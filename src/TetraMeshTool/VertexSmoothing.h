#ifndef TETRA_MESH_TOOL_VERTEX_SMOOTHING_
#define TETRA_MESH_TOOL_VERTEX_SMOOTHING_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/top.h>
#include <TetraMeshTool/QualityCalculator.h>
#include <TetraMeshTool/Quadric.h>
#include <TetraMeshTool/Journal.h>
#include <map>

namespace VolumeMesh
{

	/* store mappings from tags to vertex types */
	struct SmoothVertex
	{
		PointHandle handle;       /* point index in a tetramesh*/
		int kind;                /* the kind of vector this is (FREEVERTEX, FACETVERTEX, etc) */
		TetraMesh::Normal vec;   /* a vector associated with the vertex. for FACETVERTEX,
						               this is the normal to the plane that the vertex can move in.
						               for SEGMENTVERTEX, this is a vector in the direction of the
						               segment. */
		double val;

		bool operator <(const SmoothVertex &sv) const
		{
			return val < sv.val;
		}
		
		bool operator >(const SmoothVertex &sv) const
		{
			return val > sv.val;
		}
	};

	/* structure for holding quality and gradient information for non-smooth
	optimization-based vertex smoothing */
	struct OptTet
	{
		OptTet()
		{
			handle = TetraHandle(-1);
			volume = 0;
			volumegrad = TetraMesh::Normal(0,0,0);
			for (int i = 0; i < 6; i++)
			{
				sine[i] = 0;
				sinegrad[i] = TetraMesh::Normal(0,0,0);
			}
			rnrr = 0;
			rnrrgrad = TetraMesh::Normal(0,0,0);
			vlrms3r = 0;
			vlrms3rgrad = TetraMesh::Normal(0,0,0);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					edgelength[i][j] = 0;
					edgegrad[i][j] = TetraMesh::Normal(0,0,0);
				}
			}
			for (int i = 0; i < 4; i ++)
			{
				facearea[i] = 0;
				facegrad[i] = TetraMesh::Normal(0,0,0);
			}
		}

		TetraHandle handle;                 /* tetra indx*/
		double volume;                      /* volume of tetrahedron */
		TetraMesh::Normal volumegrad;       /* the gradient of the volume of the tet wrt vtx1 */
		double sine[6];                     /* sine of each dihedral angle of the tet */
		TetraMesh::Normal sinegrad[6];      /* gradient of each sine */
		double rnrr;                        /* root normalized radius ratio */
		TetraMesh::Normal rnrrgrad;         /* gradient of root normalized radius ratio */
		double vlrms3r;                     /* V / lrms^3 ratio */
		TetraMesh::Normal vlrms3rgrad;      /* gradient thereof */
		double edgelength[3][4];            /* the lengths of each of the edges of the tet */
		TetraMesh::Normal edgegrad[3][4];   /* the gradient of each edge length wrt vtx1 */
		double facearea[4];                 /* areas of the faces of the tet */
		TetraMesh::Normal facegrad[4];      /* the gradient of each of the face areas wrt vtx1 */
	};

	class VertexSmoother
	{
	public:
		VertexSmoother(){}
		VertexSmoother(TetraMesh *tmesh_, int qualitymetric_):tmesh(tmesh_), qualitymetric(qualitymetric_), 
			                                                  quadricVec(new QuadricContainer(tmesh_, qualitymetric_)){}
		~VertexSmoother(){}

	public:
		void setMesh(TetraMesh *tmesh_)
		{
			tmesh = tmesh_;
		}

		void setStats(ImproveStats *stats_)
		{
			stats = stats_;
		}

		void setImproveBehavior(ImproveBehavior *improvebehave_)
		{
			improvebehave = improvebehave_;
		}

		void setQualityMetric(int qualitymetric_)
		{
			qualitymetric = qualitymetric_;
		}

		/*set tools*/
		void setJournal(Journal *journals_)
		{
			journals = journals_;
		}

		void setQualityCalculator(QualityCalculator *qualitycal_)
		{
			qualitycal = qualitycal_;
		}

		void setQuadricContainer(QuadricContainer *quadricVec_)
		{
			quadricVec = quadricVec_;
		}
		
		std::map<TetraMesh::HalfFaceHandle, int> faceSegmentMap()
		{
			return facesegmap;
		}

		/* initialize the smoother */
		void initialize(int do_nonsmooth_ = 1, int fixedsmooth_ = 1);

		/* optimization-based smoothing for a single vertex */
		bool nonsmoothsinglevertex(PointHandle ph, double &worstout, int kinds);

		/* perform non-smooth optimization of a vertex's position.
		*  ph - the point to be smoothed.
		*  If we end up moving the point, return true.
		*  If for some reason we can't, return false.
		*/
		bool nonsmooth(PointHandle ph, std::vector<OptTet> &incidenttets, double &outworstqual, int smoothkinds);

		/* get vertex information: handle, vertexkind, related vec
		*/
		void vertexinformation(SmoothVertex &vinfo, PointHandle p);

		/* given two values a and b and their gradients, compute the 
		gradient of their product grad(a*b) */
		TetraMesh::Normal gradproduct(double a, double b, TetraMesh::Normal grada, TetraMesh::Normal gradb)
		{
			return TetraMesh::Normal(b*grada[0]+a*gradb[0], b*grada[1]+a*gradb[1], b*grada[2]+a*gradb[2]);
		}

		/* given two values top and bottom and their gradients, compute the 
		gradient of their quotient grad(top / bottom) */
		TetraMesh::Normal gradquotient(double top, double bot, TetraMesh::Normal gradtop, TetraMesh::Normal gradbot)
		{
			double denom = bot * bot;
			TetraMesh::Normal tmp = TetraMesh::Normal(bot*gradtop[0]-top*gradbot[0], bot*gradtop[1]-top*gradbot[1], bot*gradtop[2]-top*gradbot[2]);
			return TetraMesh::Normal(tmp[0]/denom, tmp[1]/denom, tmp[2]/denom);
		}

		/* project the vector u onto the vector v */
		TetraMesh::Normal vproject(TetraMesh::Normal u, TetraMesh::Normal v)
		{
			v.normalize();
			double tmp = u | v;
			if (v.norm() > 0)
				return TetraMesh::Normal(tmp*v[0], tmp*v[1], tmp*v[2]);
			return TetraMesh::Normal(0.0, 0.0, 0.0);
		}

		/* project the vector u onto the plane through the origin with normal v */
		TetraMesh::Normal vprojecttoplane(TetraMesh::Normal u, TetraMesh::Normal v)
		{
			return u - vproject(u,v);
		}

		/* Do the segmentation of the model
		 * Mark the half faces with the segment they are in
		 */
		void segment();


		/* Check if a half edge is sharp
		 */
		bool is_sharpedge(TetraMesh::HalfEdgeHandle he1, TetraMesh::HalfEdgeHandle he2);

		/* find the mate half edge on the boundary */
		TetraMesh::HalfEdgeHandle boudary_mate_half_edge_handle(TetraMesh::HalfEdgeHandle he);

		/* get the information about this tet needed for non-smooth
		optimization of the current quality measure */
		void getoptinfo(OptTet *opttet, SmoothVertex vinfo);

		/* compute Z, a quantity associated with circumradius computation
		TODO this code is lifted from Jonathan's tetcircumcenter computation
		in primitives.c */
		double getZ(TetraMesh::Normal tetorg, TetraMesh::Normal tetdest,
			                        TetraMesh::Normal tetfapex, TetraMesh::Normal tettapex);

		/* given a set of tets incident to a vertex, and the quality
		of the worst quality function that varies with that vertex,
		compute the active set A of quality functions very near
		the worst */
		void getactiveset(SmoothVertex sv, std::vector<OptTet> &incidenttets, 
			              std::vector<TetraMesh::Normal> &activegrads, double worstqual);


		/* finds the point on the convex hull of P nearest the origin */
		void minconvexhullpoint(PointHandle ph, std::vector<TetraMesh::Normal> P, TetraMesh::Normal &nearest);


		/* returns the basis B of S union M. S is a set of points known
		to be in the basis */
		void findbasis(std::vector<TetraMesh::Normal> M,  std::vector<Point> &S, std::vector<Point> &B);


		/* for our initial step size, we use the distance
		to the next intersection with another quality funtion.
		this is the point at which the other quality function
		becomes the worst. we use a single-term taylor expansion
		to approximate all of the quality functions as lines,
		so we'll have to do a line search to find our ultimate
		step size. */
		double getinitialalpha(std::vector<OptTet> &incidenttets, SmoothVertex vinfo, TetraMesh::Normal d, double r, double worstqual);


		/* find the best step to take to improve all of the quality
		functions affected by a vertex vtx in the search direction
		d with an expected rate of improvement r */
		void nonsmoothlinesearch(SmoothVertex vinfo, TetraMesh::Normal d, double inworstqual, double & alpha,
                                 double r, std::vector<OptTet> &incidenttets);


		/* given a boundary vertex and a list of all of it's
		incident tets, determine if it lies in an input facet */
		bool facetvertex(SmoothVertex &vinfo);

		/* given a boundary vertex and a list of all of it's
		incident tets, determine if it lies in an input segment */
		bool segmentvertex(SmoothVertex &vinfo);

		/* perform a pass of combination Laplacian / optimization-based smoothing. */
		bool smoothpass(std::vector<PointHandle> *points, double threshold, double bestmeans[], 
			            double meanqualafter[], double *minqualafter, int smoothkinds);

		/* combination lap/opt smoothing for a single point ph */
		bool smoothsinglevertex(PointHandle ph, double threshold, double &worstout, int kinds,
			                    int &lapattempts, int &lapsuccesses, int &optattempts, int &optsuccesses/*, struct arraypoolstack *outstack*/);

		/* perform a pass of Laplacian smoothing. */
		bool lapsmooth(PointHandle ph, double threshold, double &worstout, int kinds);
	private:
		TetraMesh *tmesh;
		int regionnum;                                         /*segment counter*/
		std::map<TetraMesh::HalfFaceHandle, int> facesegmap;   /*face segment map*/
		int qualitymetric;
		double sinewarpfactor;

		/* tools */
		QualityCalculator *qualitycal;
		QuadricContainer *quadricVec;
		Journal *journals;
		
		/* Journal data */
		std::vector<Point> pointvec;
		std::vector<TetraHandle> thvec;
		std::vector<PointHandle> phvec;
		std::vector<VertexHandle> vhvec;

		/*smoothing data*/
		double quadricoffset;
		double quadricscale;

		/*smoothing options*/
		int do_nonsmooth;
		int fixedsmooth;
		int usequadrics;

		/* global improvement behavior struct */
		ImproveBehavior *improvebehave;

		/* global statistics */
		ImproveStats *stats;
	};
}
#endif