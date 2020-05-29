#ifndef _ADAPTIVE_FLIP_
#define _ADAPTIVE_FLIP_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/QualityCalculator.h>

namespace VolumeMesh
{
	static const double EPS = 1e-10;
	static const double INF = 1e9;
	//static const double PI = atan2(0.0, -1.0);

	class AdaptiveFlip
	{
	public:
		AdaptiveFlip():mesh(NULL){}
		AdaptiveFlip(TetraMesh *mesh_):mesh(mesh_){qualcalculator.setMesh(mesh_);}
		~AdaptiveFlip(){}

		void setMesh(TetraMesh *mesh_)
		{
			mesh = mesh_;
			qualcalculator.setMesh(mesh_);
		}

		// situation observation
		void situationObservation(int s[2], std::vector<TetraHandle> &tetvec);

		double vertexTriangleDistance(Point triangle[3], Point p, int & pStatus);
		int AdaptiveFlip::checkPointInsideTriangle(Point triangle[3], Point p);
		void findcrossedges(TetraMesh::HedronHandle th_, TetraMesh::HalfEdgeHandle &e1, TetraMesh::HalfEdgeHandle &e2);
		bool kitesituationcheck(TetraHandle th_, TetraMesh::HalfEdgeHandle e1, TetraMesh::HalfEdgeHandle e2);
		bool trianglesituationcheck(TetraHandle th_, TetraMesh::VertexHandle vh_);

		void kiteAdaptiveFlip(TetraHandle th_, TetraMesh::HalfEdgeHandle e1, TetraMesh::HalfEdgeHandle e2);
		void triangleAdaptiveFlip(TetraHandle th_, TetraMesh::VertexHandle vh_);

	private:
		TetraMesh *mesh;
		QualityCalculator qualcalculator;
	};

}
#endif