#ifndef _TETRAMESH_TOOL_QUALITY_CALCULATOR
#define _TETRAMESH_TOOL_QUALITY_CALCULATOR

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/top.h>

namespace VolumeMesh
{
	class QualityCalculator
	{
	public:
		QualityCalculator(){}
		QualityCalculator(TetraMesh *tmesh_):tmesh(tmesh_){}
		~QualityCalculator(){}

	public:
		double getZ(Point tetorg, Point tetdest, Point tetfapex, Point tettapex);
		double minsine(Point p1, Point p2, Point p3, Point p4);
		double biasedminsine(Point p1, Point p2, Point p3, Point p4);
		double meansine(Point p1, Point p2, Point p3, Point p4);
		double minsineandedgeratio(Point p1, Point p2, Point p3, Point p4);
		double radiusratio(Point p1, Point p2, Point p3, Point p4);
		double vlrms3ratio(Point p1, Point p2, Point p3, Point p4);
		double warpedminsine(Point p1, Point p2, Point p3, Point p4);
		double minmaxangle(Point p1, Point p2, Point p3, Point p4, bool max);
		double tetquality(Point p1, Point p2, Point p3, Point p4, int qualityMetric);
		double minstackquality(TetraMesh *mesh, const std::vector<TetraHandle> &tetstack, int qualmeasure);
		double minstackquality(TetraMesh *mesh, TetraStack &tetstack, int qualmeasure);
		double minmeshquality(TetraMesh *mesh, int qualmetric);
		void meshquality(TetraMesh *mesh, int qualmeasure, double *minqual, TetraStack *tetstac=NULL);
		void fillstackqual(TetraMesh *mesh, int qualmeasure, double *minqual, double threshold=0, TetraStack *tetstack=NULL);
		void stackquality(TetraMesh *mesh, TetraStack &tetstack, int qualmeasure, double meanqual[], double *minqual);
		void stackquality(TetraMesh *mesh, std::vector<TetraHandle> &tetstack, int qualmeasure, double meanqual[], double *minqual);
		double worstinputangle(TetraMesh *mesh);
		double worstquality(TetraMesh *mesh);
		double worstquality(const std::vector<TetraHandle> &tetvec);
		double getboundaryedgeangle(TetraMesh *mesh, TetraHandle th, int v1, int v2, int vleft, int vright);

		void setMesh(TetraMesh *mesh_)
		{
			tmesh = mesh_;
		}

		void setQualityMetric(int qualmetric_)
		{
			qualmetric = qualmetric_;
		}

		void setSineWarpFactor(double sinewarpfactor_)
		{
			sinewarpfactor = sinewarpfactor_;
		}

		/* convert from a sin of an angle to that angle in degrees (does not test obtuseness) */
		double sintodeg(double insine)
		{
			return asin(insine) * 180.0 / PI;
		}

		/* convert from an angle in degrees to sin */
		double degtosin(double inangle)
		{
			return sin((inangle) * (PI / 180.0));
		}

		double tantodeg(double intan)
		{
			return atan(intan) * 180.0 / PI;
		}

		double radtodeg(double inangle)
		{
			return (inangle * 180) / PI;
		}

	private:
		TetraMesh *tmesh;
		int qualmetric;
		double sinewarpfactor;
	};
}
#endif