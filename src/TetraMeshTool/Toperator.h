#ifndef _VOLUME_MESH_TOPOGLOGYICAL_OPERATOR_
#define _VOLUME_MESH_TOPOGLOGYICAL_OPERATOR_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <VolumeMesh/Mesh/HexadMesh.h>
#include <TopologyProcess/TetraHedronEntity.h>
#include <TetraMeshTool/top.h>
#include <TetraMeshTool/Journal.h>
#include <TetraMeshTool/QualityCalculator.h>
#include <TetraMeshTool/StarBase.h>

namespace VolumeMesh
{
	class TMeshToperator
	{
	public:
		TMeshToperator():edgeremoval(1),singlefaceremoval(1), multifaceremoval(1)
		{
		}
		TMeshToperator(TetraMesh *mesh_)
		{
			tmesh = mesh_;
		}
		~TMeshToperator(){}

		void setMesh(TetraMesh *mesh_)
		{
			tmesh = mesh_;
		}

		void setQualityMetric(int qualityMetric_)
		{
			qualityMetric = qualityMetric_;
		}

		void setJournal(Journal *journal_)
		{
			optjournal = journal_;
		}

		void setStats(ImproveStats *stats_)
		{
			stats = stats_;
		}

		void setImproveBehavior(ImproveBehavior *improvebehave_)
		{
			improvebehave = improvebehave_;
		}

		void setQualityCalculator(QualityCalculator *qualCalculator_)
		{
			qualcalculator = qualCalculator_;
		}

		void getRecoverData(int type, std::vector<Point> &pointvec, std::vector<TetraHandle> &thvec, 
			std::vector<PointHandle> &phvec, std::vector<VertexHandle> &vhvec);

		void setRecoverData(int type, std::vector<Point> &pointvec, std::vector<TetraHandle> &thvec, 
			                std::vector<PointHandle> &phvec, std::vector<VertexHandle> &vhvec);

		bool recover(int optType, std::vector<TetraHandle> *tetraVec_ = NULL);


		/* attempt to remove the edge he_. do flip22 or flip32.
		   [PARAM] he_ : is the edge to be removed
		   [PARAM] oldminqual : original worst quality
		   [PARAM] newtets : new generated tetras
		   [PARAM] outminqual : the worst quality after edge removal
		   [PARAM] boundary : indicate if the edge is a boundary one
		   [RETRUN] bool : TRUE edge removal success
		                   FALSE failed
		   */
		bool removeedge(HalfEdgeHandle he_, double &oldminqual, std::vector<TetraHandle> &newtets, 
			            double &outminqual, bool &boundary);

		/* attempt to remove the face hf_. do flip23.
		   [PARAM] hf_ : is the face to be removed
		   [PARAM] oldminqual : original worst quality
		   [PARAM] newtets : new generated tetras
		   [PARAM] outminqual : the worst quality after edge removal
		   [PARAM] boundary : indicate if the edge is a boundary one
		   [RETRUN] bool : TRUE edge removal success
		                   FALSE failed
		   */
		bool removeface(HalfFaceHandle hf_, double oldminqual, std::vector<TetraHandle> &newtets,
			            double &outminqual, bool &boundary);


		/* fill Q and K tables for Klincsek's algorithm */
		void filltables(HalfEdgeHandle he_, int *ring, int ringcount, double oldminqual, 
			            double Q[][MAXRINGTETS], int K[][MAXRINGTETS]);

		/* perform a pass of topological improvement.
		for now, this means trying to remove each edge
		of each tet in the stack that is passed in,
		and if no edge can be removed, trying to remove
		each face. */
		bool topologyPass(std::vector<TetraHandle> *tetstack, std::vector<TetraHandle> *outstack, double &minqualafter);
		
		/* go after the worst tets with contraction */
		void contractworst(double percentinsert, double bestmeans[], double outmeanqual[], double *outminqual, bool desperate);

		/* for each tet in the stack, try to contract its edges */
		bool contractpass(std::vector<TetraHandle> *tetstack, std::vector<TetraHandle> *outstack, 
			              double &minqualafter, bool justfirstedge);
		/* try edge removal on all 6 edges of the specified tet. if removal
		succeeds, return 1. Otherwise, return 0. */
		bool tryedgeremove(TetraHandle tetrahandle, std::vector<TetraHandle> &newtetras, double &minafterqual);

		/* try a 2-3 flip on all four faces of the specified tet
		returns 1 if a flip was performed, 0 otherwise */
		bool tryremovefaces(TetraHandle tetrahandle, std::vector<TetraHandle> &newtetras, double &minafterqual);

		/* try to contract all six edges of a tet */
		bool tryedgecontract(TetraHandle tetrahandle, std::vector<TetraHandle> *outstack, double *minqualafter);
		/* given a mesh and a percentage p,
		return the worst numtets * p tets in the mesh */
		void fillstackpercent(TetraMesh *mesh, std::vector<TetraHandle> *stack, int qualmeasure,
			                  double percent);

		void fillstackqual(TetraMesh *mesh, std::vector<TetraHandle> *stack, int qualmeasure,
			               double threshold);

	public:
		/************************************************************************/
		/* Topology operation                                                   */
		/************************************************************************/
		int flip23(HalfFaceHandle &hf_, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip32(HalfEdgeHandle &he_, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip22(HalfEdgeHandle &he_, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip13(HalfFaceHandle &hf_, Point &insertPoint, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip31(HalfFaceHandle hf[3], PointHandle &midPoint, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip14(TetraHandle &h_, Point &insertPoint, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip41(PointHandle &ph_, std::vector<TetraHandle> *tetraVec_ = NULL);
		int flip12(HalfEdgeHandle &he_, PointHandle &ph_, std::vector<TetraHandle> *tetraVec_ = NULL);
		int edge_contract(HalfEdgeHandle &he_, PointHandle *newpoint = NULL);

	public:
		/************************************************************************/
		/* Recover operation                                                    */
		/************************************************************************/
		int recover_flip12(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip13(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip23(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip22(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip32(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip31(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip14(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_flip41(std::vector<TetraHandle> *tetraVec_ = NULL);
		int recover_edge_contract(std::vector<TetraHandle> *tetraVec_ = NULL);

		PointHandle testingpoint()
		{
			return ec_contract_point;
		}

	public:
		//-------------------------------------control parameter-------------------------------------//
		int edgeremoval;
		int singlefaceremoval;
		int multifaceremoval;

	private:
		TetraMesh *tmesh;
		int qualityMetric;
		Journal *optjournal;
		QualityCalculator *qualcalculator;

		/* global improvement behavior struct */
		ImproveBehavior *improvebehave;

		/* global statistics */
		ImproveStats *stats;

		/* Journal data */
		std::vector<Point> jpointvec;
		std::vector<TetraHandle> jthvec;
		std::vector<PointHandle> jphvec;
		std::vector<VertexHandle> jvhvec;

		//-------------------------------------flip2-3 recover data------------------------------------//
		bool rec_flip23_valid;
		PointHandle rec_flip23_top;          /**< top point of the five points  */
		PointHandle rec_flip23_bot;          /**< bottom point of the five points  */
		PointHandle rec_flip23_p[3];         /**< the other three points  */
		TetraHandle rec_flip23_old_tetra[2]; /**< the two original tetrahedrons */
		TetraHandle rec_flip23_tetra[3];     /**< the three new tetrahedrons  */

		//-------------------------------------flip2-2 recover data------------------------------------//
		bool rec_flip22_valid;
		PointHandle rec_flip22_top;          /**< the point is not in either of the coplanar halffaces  */
		PointHandle rec_flip22_p1;           /**< endpoint of the input halfedge  */
		PointHandle rec_flip22_p2;           /**< endpoint of the input halfedge  */
		PointHandle rec_flip22_p[2];         /**< the other two points  */
		TetraHandle rec_flip22_old_tetra[2]; /**< the two original tetrahedorns */
		TetraHandle rec_flip22_tetra[2];     /**< the two new tetrahedrons  */

		//-------------------------------------flip3-2 recover data------------------------------------//
		bool rec_flip32_valid;
		PointHandle rec_flip32_top;          /**< top point of the five points  */
		PointHandle rec_flip32_bot;          /**< bottom point of the five points  */
		PointHandle rec_flip32_p[3];         /**< the other three points  */
		TetraHandle rec_flip32_old_tetra[3]; /**< the three original tetrahedrons */
		TetraHandle rec_flip32_tetra[2];     /**< the two new tetrahedrons  */

		//-------------------------------------flip1-3 recover data------------------------------------//
		bool rec_flip13_valid;
		PointHandle rec_flip13_p[4];         /**< the point handles of the original tetra  */
		TetraHandle rec_flip13_old_tetra;    /**< the original tetrehedron */
		TetraHandle rec_flip13_tetra[3];     /**< the three new tetrahedrons  */

		//-------------------------------------flip3-1 recover data------------------------------------//
		bool rec_flip31_valid;
		Point rec_flip31_midpoint;           /**< the midpoint shared by three boundary faces */
		PointHandle rec_flip31_top;          /**< top point of the five points  */
		PointHandle rec_flip31_face_p[3];    /**< the three points of divided face  */
		TetraHandle rec_flip31_old_tetra[3]; /**< the three original tetrahedrons */
		TetraHandle rec_flip31_tetra;        /**< the new tetrahedron  */

		//-------------------------------------flip1-4 recover data------------------------------------//
		bool rec_flip14_valid;
		PointHandle rec_flip14_p[4];        /**< the point handles of the original tetra  */
		TetraHandle rec_flip14_old_tetra;   /**< the original tetrahedron */
		TetraHandle rec_flip14_tetra[4];    /**< the four new tetrahedrons  */

		//-------------------------------------flip4-1 recover data------------------------------------//
		bool rec_flip41_valid;
		Point rec_flip41_midpoint;           /**< the midpoint shared by three boundary faces */
		TetraHandle rec_flip41_old_tetra[4]; /**< the four original tetrahedrons */
		TetraHandle rec_flip41_tetra;        /**< the new tetrahedron */

		//-------------------------------------flip1-2 recover data------------------------------------//
		bool rec_flip12_valid;
		PointHandle rec_flip12_p[4];         /**< the point handles of the original tetra  */
		TetraHandle rec_flip12_old_tetra;    /**< the original tetrahedron */
		TetraHandle rec_flip12_tetra[2];     /**< the two new tetrahedrons  */

		//-----------------------------------edge contract recover data--------------------------------//
		bool rec_ec_valid;
		Point rec_ec_erase_point;            /**< coordinate of erased point */
		PointHandle rec_ec_point_handle;     /**< the point handle of the reserved endpoint of contracted edge */
		Point rec_ec_point;                  /**< the coordinate of the reserved endpoint of contracted edge */
		std::vector<PointHandle> rec_ec_tetravertex;   /**< the vector contains the other points (except the two endpoints of the edge) of each erased tetras */
		std::vector<TetraHandle> rec_ec_old_tetra;     /**< the original tetrahedrons incident to the contracting edge */
		std::vector<TetraHandle> rec_ec_to_point_vertex_container;    /**< set of VertexHandles that their original PointHanle is to_point_handle  */

		//TESTING 
		PointHandle ec_contract_point;
	};
}
#endif