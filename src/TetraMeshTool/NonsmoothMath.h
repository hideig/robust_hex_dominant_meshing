#ifndef NONSMOOTH_MATH_
#define  NONSMOOTH_MATH_
#include <TetraMeshTool/NonsmoothTypes.h>
#include < math.h>
namespace VolumeMesh
{
	inline NonsmoothPoint gradproduct(double a, double b, NonsmoothPoint grada, NonsmoothPoint gradb)
	{
		return a * gradb + b * grada;
	}
	inline NonsmoothPoint gradquotient(double top, double bot, NonsmoothPoint gradtop, NonsmoothPoint gradbot)
	{
		double denom = bot * bot;
		NonsmoothPoint tmp = bot *gradtop - top * gradbot;
		return tmp / denom;
	}
	// count dot product for two vector
	inline double dotProduct(double _point_1[3], double _point_2[3])
	{
		return _point_1[0] * _point_2[0] + _point_1[1] * _point_2[1] + _point_1[2] * _point_2[2];
	}
	//count cross product for two vector
	inline void crossProduct(double _point_1[3], double _point_2[3], double result[3])
	{
		result[0] = _point_1[1]*_point_2[2] - _point_1[2]*_point_2[1];
		result[1] = _point_1[2]*_point_2[0] - _point_1[0]*_point_2[2];
		result[2] = _point_1[0]*_point_2[1] - _point_1[1]*_point_2[0];
	}
	//count the length of a vector
	//double lengths (NonsmoothPoint _p)
	//{
	//	//return sqrt(pow(_p.point_coord[0], 2.0) + pow(_p.point_coord[1], 2.0) + pow(_p.point_coord[2], 2.0));
	//	return sqrt(_p |_p);
	//}

	//count the distance between two point
	inline double distance (NonsmoothPoint _p_1, NonsmoothPoint _p_2)
	{
		NonsmoothPoint p;
		p = _p_1 - _p_2;
		return sqrt(p | p);
	}
	//get the unit vector of vector
	inline NonsmoothPoint unitVector(NonsmoothPoint _v)
	{
		return _v / sqrt(_v | _v);
	}
	// project the vector v_1 to the vector v_2
	inline NonsmoothPoint project (NonsmoothPoint _v_1, NonsmoothPoint _v_2)
	{
		if (sqrt(_v_2 | _v_2) > 0)
		{
			double scale = (_v_1 | _v_2) / sqrt(_v_2 | _v_2);
			NonsmoothPoint direction = unitVector(_v_2);
			return scale * direction;
		}
		else
			return 0 * _v_2;
	}
	// project the vector v to the plane with the normal n
	inline NonsmoothPoint projectToPlane (NonsmoothPoint _v, NonsmoothPoint _n)
	{
		return _v - project(_v, _n);
	}
	//count the volume of the tetrahedron
	inline double TetrahedronVolume(OneTetrahedron _tetrahedron)
	{
		NonsmoothPoint p[4];
		p[0] = _tetrahedron.main_point;
		p[1] = _tetrahedron.minor_points[0];
		p[2] = _tetrahedron.minor_points[1];
		p[3] = _tetrahedron.minor_points[2];
		return abs((((p[1]-p[0]) % (p[2]-p[0])) | (p[3]-p[0])) / 6);
	}
	inline double TetrahedronVolume(NonsmoothPoint p[4])
	{
		return abs((((p[1]-p[0]) % (p[2]-p[0])) | (p[3]-p[0])) / 6);
	}
	/* check if a facet (indicated by fp1, fp2, fp3) is oriented to a vertex v */
	inline int orient(NonsmoothPoint v, NonsmoothPoint fp1, NonsmoothPoint fp2, NonsmoothPoint fp3)
	{
		NonsmoothPoint center;
		NonsmoothPoint facenormal;
		double result;
		center = (fp1 + fp2 + fp3) / 3.0;
		facenormal = (fp2 - fp1) % (fp3 - fp1);
		facenormal = unitVector(facenormal);
		result = (v - center) | facenormal;
		if (result > 0)
			return 1;
		if (result == 0)
			return 0;
		return -1;
	}
	inline double SinToDegree(double _value)
	{
		return asin(_value) * 180 / PI;
	}
}
#endif