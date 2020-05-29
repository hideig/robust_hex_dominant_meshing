#ifndef TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_QUALITY_
#define TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_QUALITY_

#include <TetraMeshTool/NonsmoothTypes.h>


namespace VolumeMesh
{
	class NonsmoothQuality
	{
	public:
		NonsmoothQuality(){}
		NonsmoothQuality(int _quality_kind):quality_kind_(_quality_kind){}
		~NonsmoothQuality(){}
		void setTetrahedron(OneTetrahedron _tetrahedron)
		{
			tetrahedron_ = _tetrahedron;
		}
		void setTetrahedron2(OneTetrahedron2 _tetrahedron2)
		{
			tetrahedron2_ = _tetrahedron2;
			tetrahedron_.main_point_index = tetrahedron2_.four_points_index[0];
			tetrahedron_.minor_points_index[0] = tetrahedron2_.four_points_index[1];
			tetrahedron_.minor_points_index[1] = tetrahedron2_.four_points_index[2];
			tetrahedron_.minor_points_index[2] = tetrahedron2_.four_points_index[3];
		}
		void setQualityKind(int _quality_kind)
		{
			quality_kind_ = _quality_kind;
		}
		double getQuality(NonsmoothPoint p[4],int _quality_kind)
		{
			NonsmoothPoint p1,p2,p3,p4;
			p1 = p[0];
			p2 = p[1];
			p3 = p[2];
			p4 = p[3];
			switch(_quality_kind)
			{
			case MINSINE: return getMinsine(p);break;
			case VOLLENGTH: return getVomlength(p);break;
			case RADIUSRATIO: return getRadiusratio(p);break;
			case MINSINE2:return getMinsine2(p);break;
			case BIASEDMINSINE:return getBiasedMinsine(p);break;
			default:
				return 0;
			}
		}
		double getMinsine(NonsmoothPoint p[4]);
		double getVomlength(NonsmoothPoint p[4]);
		double getRadiusratio(NonsmoothPoint p[4]);
		double getMinsine2(NonsmoothPoint p[4]);
		double getBiasedMinsine(NonsmoothPoint p[4]);
		double getZ(NonsmoothPoint p[4]);
	protected:
	private:
		int quality_kind_;
		OneTetrahedron tetrahedron_;
		OneTetrahedron2 tetrahedron2_;
	};
}
#endif