#ifndef TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_TYPE_
#define TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_TYPE_
#include "ConstantData.h"
#include <set>
namespace VolumeMesh
{

	//define point structure
	struct NonsmoothPoint
	{
		double point_coord[3];
		bool boundary;
		NonsmoothPoint(){};
		NonsmoothPoint(double p[3])
		{
			point_coord[0] = p[0];
			point_coord[1] = p[1];
			point_coord[2] = p[2];
		}
		NonsmoothPoint(double p1, double p2, double p3)
		{
			point_coord[0] = p1;
			point_coord[1] = p2;
			point_coord[2] = p3;
		}
		friend NonsmoothPoint operator + (const NonsmoothPoint & left, const NonsmoothPoint & right)
		{
			NonsmoothPoint result;
			result.point_coord[0] = left.point_coord[0] + right.point_coord[0];
			result.point_coord[1] = left.point_coord[1] + right.point_coord[1];
			result.point_coord[2] = left.point_coord[2] + right.point_coord[2];
			return result;
		}
		friend NonsmoothPoint operator - (const NonsmoothPoint & left, const NonsmoothPoint & right)
		{
			NonsmoothPoint result;
			result.point_coord[0] = left.point_coord[0] - right.point_coord[0];
			result.point_coord[1] = left.point_coord[1] - right.point_coord[1];
			result.point_coord[2] = left.point_coord[2] - right.point_coord[2];
			return result;
		}
		friend double operator | (const NonsmoothPoint & left, const NonsmoothPoint & right)
		{
			return left.point_coord[0] * right.point_coord[0] + left.point_coord[1] * right.point_coord[1] + left.point_coord[2] * right.point_coord[2];

		}

		friend NonsmoothPoint operator % (const	NonsmoothPoint	 & left, const NonsmoothPoint & right)
		{
			NonsmoothPoint new_point;
			new_point.point_coord[0] = left.point_coord[1] * right.point_coord[2] - left.point_coord[2] * right.point_coord[1];
			new_point.point_coord[1] = left.point_coord[2] * right.point_coord[0] - left.point_coord[0] * right.point_coord[2];
			new_point.point_coord[2] = left.point_coord[0] * right.point_coord[1] - left.point_coord[1] * right.point_coord[0];
			return new_point;
		}
		friend NonsmoothPoint operator * (const NonsmoothPoint & left, const double & right)
		{
			NonsmoothPoint new_point;
			new_point.point_coord[0] = left.point_coord[0] * right;
			new_point.point_coord[1] = left.point_coord[1] * right;
			new_point.point_coord[2] = left.point_coord[2] * right;
			return new_point;
		}
		friend NonsmoothPoint operator * ( const double & left, const NonsmoothPoint & right)
		{
			NonsmoothPoint new_point;
			new_point.point_coord[0] = right.point_coord[0] * left;
			new_point.point_coord[1] = right.point_coord[1] * left;
			new_point.point_coord[2] = right.point_coord[2] * left;
			return new_point;
		}
		friend NonsmoothPoint operator / (const NonsmoothPoint & top, const double & bottom)
		{
			return top * (1 / bottom);
		}
		friend bool operator == (const NonsmoothPoint & left, const NonsmoothPoint & right)
		{
			return (( left.point_coord[0] == right.point_coord[0] ) && ( left.point_coord[1] == right.point_coord[1] ) && ( left.point_coord[2] == right.point_coord[2] ));
		}
		friend bool operator < (const NonsmoothPoint & left, const NonsmoothPoint & right)
		{
			return ((left.point_coord[0] < right.point_coord[0]) || ((left.point_coord[0] == right.point_coord[0]) &&(left.point_coord[1] < right.point_coord[1])) ||
				((left.point_coord[0] == right.point_coord[0]) && (left.point_coord[1] == right.point_coord[1]) && (left.point_coord[2] < right.point_coord[2])));
		}
		friend NonsmoothPoint operator -(NonsmoothPoint & left)
		{
			left.point_coord[0] = - left.point_coord[0];
			left.point_coord[1] = - left.point_coord[1];
			left.point_coord[2] = - left.point_coord[2];
			return left;
		}
		double operator [](int index)
		{
			return point_coord[index];
		}
		/*friend std::ostream &operator << (std::ostream & os, NonsmoothPoint & p)
		{
			os<<p.point_coord[0]<<" "<<p.point_coord[1]<<" "<<p.point_coord[2];
			return os;
		}*/
	};
	struct OneTetrahedron
	{
		NonsmoothPoint main_point;
		int main_point_index;
		NonsmoothPoint minor_points[3];
		int minor_points_index[3];
	};
	struct OneTetrahedron2
	{
		//main_point_index used for tagging the vertices
		//int main_point_index;
		int four_points_index[4];
		bool boundary;
	/*	OneTetrahedron2 operator = (OneTetrahedron2 & right)
		{
			four_points_index[0] = right.four_points_index[0];
			four_points_index[1] = right.four_points_index[1];
			four_points_index[2] = right.four_points_index[2];
			four_points_index[3] = right.four_points_index[3];
		}*/
	};
	struct IncidentPoints
	{
		int incident_point_num;
		int incident_points[MAXINCIDENTPOINT];
	};
	struct IncidentTetrahedrons
	{
		int incident_tetrahedron_num;
		OneTetrahedron incident_tetrahedrons[MAXINCIDENTHETRDRON];
	};
	struct HedronGradient
	{
		HedronGradient()
		{
			volume = 0;
			volumegrad = NonsmoothPoint(0,0,0);
			for (int i = 0; i < 6; i++)
			{
				sine[i] = 0;
				sinegrad[i] = NonsmoothPoint(0,0,0);
			}
		}

		double volume;                      /* volume of tetrahedron */
		NonsmoothPoint volumegrad;       /* the gradient of the volume of the tet wrt vtx1 */
		double sine[6];                     /* sine of each dihedral angle of the tet */
		NonsmoothPoint sinegrad[6];      /* gradient of each sine */
	};
	struct IncidentTetrahedrons2
	{
		int incident_hedron_num;
		int incident_tetrahedrons[MAXINCIDENTHETRDRON];
		//HedronGradient hedrongradient[MAXINCIDENTHETRDRON];
	};
	
}

#endif
