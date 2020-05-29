#ifndef _TETRA_MESH_TOOL_TETRA_MESH_IMPROVER
#define _TETRA_MESH_TOOL_TETRA_MESH_IMPROVER

#include <iostream>
#include <TetraMeshTool/Journal.h>
#include <TetraMeshTool/Toperator.h>
#include <TetraMeshTool/VertexInsertion.h>
#include <TetraMeshTool/VertexSmoothing.h>
#include <TetraMeshTool/Quadric.h>
#include <TetraMeshTool/VertexNonsmooth.h>
namespace VolumeMesh
{
	class TetraMeshImprover
	{
	public:
		TetraMeshImprover():mesh(NULL){}
		TetraMeshImprover(TetraMesh *mesh_):mesh(mesh_){}
		~TetraMeshImprover(){}

		TetraMesh* getMesh()
		{
			return mesh;
		}

		void setMesh(TetraMesh *mesh_)
		{
			mesh = mesh_;
		}

		/* control parameters setting */
		void setQualityMetric(int qualityMetric_)
		{
			qualitymetric = qualityMetric_;
		}

		/* initialize the improver */
		void initialize();
		void parseimprovecommandline(ImproveBehavior *b);
		void parseimprovestateline(ImproveStats *stats);

		/* tetrahedral mesh improvement */
		bool meshImproving(TetraMesh *mesh_, double &qualafter);

		/* pre-improvement initialization code */
		void improveinit(TetraMesh *mesh, double bestmeans[NUMMEANTHRESHOLDS]);

		/* see if new means contains any better means. if so,
		update best means and return true */
		bool meanimprove(double bestmeans[], double newmeans[], int passtype);
		
		/* run a pass (smoothing, topo, insertion). return true
		if we have reached the desired quality */
		bool pass(int passtype, TetraMesh *mesh, double threshold, 
			      bool &minsuccess, bool &meansuccess,
			      int passnum, double bestmeans[]);

		void getextremeangles(TetraMesh *mesh, double &outsmallestangle, double &outbiggestangle,int stat[180]);
		void getextremeangles(TetraMesh *mesh, double &outsmallestangle, double &outbiggestangle);

		void rebuildmesh(TetraMesh *mesh);
	protected:
	private:
		TetraMesh *mesh;

		// equipments
		Journal journals;
		VertexSmoother vsmoother;
		Nonsmoother vnonsmoother;
		VertexInsert vinserter;
		TMeshToperator toperator;
		QualityCalculator qualcalculator;
		QuadricContainer quadric;

		// control parameters
		int qualitymetric;
		int edgeremoval;
		int singlefaceremoval;
		int multifaceremoval;

		// state parameters
		int journalentry;

		//improvement elements
		//std::vector<TetraHandle> improveTetraVec;
		//std::vector<PointHandle> improvePointVec;

		/* global improvement behavior struct */
		ImproveBehavior improvebehave;

		/* global statistics */
		ImproveStats stats;
	};
}
#endif