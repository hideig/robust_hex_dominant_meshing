/*****************************************************************************/
/*  TetraMesh                                                                */
/*  Vertex insertion routines                                                */
/*  Base on steller                                                          */
/*****************************************************************************/
#ifndef _TETRAMESH_VERTEX_INSETION
#define _TETRAMESH_VERTEX_INSETION

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TopologyProcess/TetraHedronEntity.h>
#include <TetraMeshTool/Toperator.h>
#include <TetraMeshTool/QualityCalculator.h>
#include <TetraMeshTool/Journal.h>
#include <TetraMeshTool/Quadric.h>
#include <TetraMeshTool/VertexSmoothing.h>
#include <TetraMeshTool/AdaptiveFlip.h>

namespace VolumeMesh
{
	/* A compare class for determining whether the tet quality is greater than the threshold*/
	template<class T1, class T2>
	class GTThreshold : public std::binary_function<T1, T2, bool>
	{
	public:
		bool operator()(const T1 &a, const T2 &threshold) const
		{
			return a.quality>threshold;
		}
	};

	static double cavitydepthtable[DEPTHTABLESIZE] = {1.0, 1.6, 2.3, 2.9, 3.3, 
		                                 3.3, 3.3, 3.3, 3.3, 3.3};
	class VertexInsert
	{
	public:
		VertexInsert()
		{
			maxcavitytet = 0;
			maxcavityface = 0;
			maxcavityedge = 0;
			maxstacktet = 0;
			maxstackface = 0;
			//maxstackedge = 0;
		}
		VertexInsert(TetraMesh *tmesh_, int qualitymetric_): tmesh(tmesh_), qualitymetric(qualitymetric_)
		{
			qualcalculator = new QualityCalculator;
			qualcalculator->setMesh(tmesh);
			maxcavitytet = 0;
			maxcavityface = 0;
			maxcavityedge = 0;
			maxstacktet = 0;
			maxstackface = 0;
			//maxstackedge = 0;
		}
		~VertexInsert(){}

	public:
		void setMesh(TetraMesh *tmesh_)
		{
			tmesh = tmesh_;
		}

		void setQualityMetric(int qualitymetric_)
		{
			qualitymetric = qualitymetric_;
		}

		// set insertion options
		void setInsertThreshold(double threshold)
		{
			insertthreshold = threshold;
		}

		void enableCavityDeleted(bool opt)
		{
			cavityconsiderdeleted = opt;
		}

		void cavityDepthLimit(int limit)
		{
			cavdepthlimit = limit;
		}

		/************************************************************************/
		/* Convenience functions                                                */
		/************************************************************************/

		/* A predicate function for cavityedge comparision */
		static bool comCavityEdge_less(const CavityEdgeorTet &e1, const CavityEdgeorTet &e2);
		/* A predicate function for cavitytet comparison*/
		static bool comCavityTet(const CavityTet &a, const CavityTet &b);
		static bool comCavityTetHandle(const CavityTet &a, const TetraHandle & handle);
		/* A predicate function for ImproveTetra comparison*/
		static bool comImproveTets(const ImproveTetra &a, const ImproveTetra &b);
		/* get four points of a tetrahedron indexed by a TetraHandle*/
		static bool tetPoints(TetraMesh *mesh, TetraHandle handle, Point *p);
		/* get three points of a half face indexed by a TetraHandle*/
		static bool halfFacePoints(TetraMesh *tmesh, TetraMesh::HalfFaceHandle hf, Point *p);
		/* get three points of a half face indexed by a TetraHandle*/
		static bool halfFacePointHandle(TetraMesh *tmesh, TetraMesh::HalfFaceHandle hf, PointHandle *p);
		void initImproveCommand();  // should be run at the beginning of improvement

	public:
		void targetTetras(double &outminqual, int *trytetcnt = NULL, int *succcnt = NULL);
		bool insertPass(TetraStack &tetstack, int qualmeasure, double &minqualafter, double okayqual, int *succcnt = NULL);
		bool insertVertex(TetraMesh::HedronHandle tethandle, std::vector<TetraMesh::HedronHandle> *outtet = NULL);
		//bool facetInsert(TetraMesh::HalfFaceHandle facehandle, PointHandle &pointhandle, 
		//	             std::vector<TetraMesh::HedronHandle> &newtets, std::vector<TetraMesh::HalfFaceHandle> &outfaces);
		bool bodyinsert(TetraMesh::HedronHandle tetrahandle, PointHandle &newp, 
			            std::vector<HalfFaceHandle> &newfaces, std::vector<TetraHandle> &newtets);
		//bool segmentinsert(TetraHandle tetrahandle, HalfEdgeHandle splitedge, PointHandle &newp, 
		//	               std::vector<HalfFaceHandle> &newfaces, std::vector<TetraHandle> &newtets);
		void buildcavitydag(PointHandle pnew, std::vector<TetraHandle> initialC, std::vector<TetraMesh::HalfFaceHandle> initialFaces, 
			                std::vector<CavityTet> &outcavity, bool allowvertexdeletion);
		unsigned int findCavityTet(std::vector<CavityTet>::iterator begin, std::vector<CavityTet>::iterator end, TetraHandle handle);
		std::vector<CavityFace>::iterator findCavityFace(std::vector<CavityFace>::iterator begin, std::vector<CavityFace>::iterator end, 
			                                             TetraMesh::HalfFaceHandle handle);
		void tetadjacencies(TetraMesh::HalfFaceHandle hf, TetraHandle *outin);
		void maxCavity(PointHandle pnew, std::vector<CavityTet> &cavity, std::vector<TetraHandle> &outtets, 
			           double *worstdeleted, double *worstincavity);
		int numParents(CavityTet &tet);
		void addCavityTetFace(CavityTet &tet, CavityFace &cface);
		void cavityLabel(std::vector<CavityTet> &cavity, unsigned int ctetidx, std::vector<int> *labeltets = NULL);
		void antiCavityLabel(std::vector<CavityTet> &cavity, unsigned int ctet_iter, std::vector<int> *labeltets = NULL);
		bool smoothInsertVertex(Point &p, std::vector<Point> &faces, int qualmetric);

		void cavityLabel_loop(std::vector<CavityTet> &cavity, unsigned int ctetidx_, std::vector<int> &labeltets);
		void antiCavityLabel_loop(std::vector<CavityTet> &cavity, unsigned int ctetidx_, std::vector<int> &labeltets);

		void getSuccedTet(std::vector<int> &vec)
		{
			vec.clear();
			copy(succvec.begin(), succvec.end(), back_inserter(vec));
		}

		void getStatisticData(int &maxcavitytet_, int &maxcavityface_, int &maxcavityedge_, int &maxstacktet_, int &maxstackface_)
		{
			maxcavitytet_ = maxcavitytet;
			maxcavityface_ = maxcavityface;
			maxcavityedge_ = maxcavityedge;
			maxstacktet_ = maxstacktet;
			maxstackface_ = maxstackface;
		}
	private:
		TetraMesh *tmesh;
		int qualitymetric;
		QualityCalculator *qualcalculator;

		/* insertion options */
		double insertthreshold;       /* percent worst tets */
		bool cavityconsiderdeleted;   /* consider enlarging cavity for deleted tets? */
		int cavdepthlimit;            /* only allow initial cavity to includes tets this deep */
		int cavdeep;

		//testing
		std::vector<int> succvec;
		int maxcavitytet;
		int maxcavityface;
		int maxcavityedge;

		int maxstacktet;
		int maxstackface;
		//int maxstackedge;
	};
	
}
#endif