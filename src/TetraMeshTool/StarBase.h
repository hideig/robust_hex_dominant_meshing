#ifndef _TOPOLOGICAL_OPERATOR_STARBASE_DEFINE
#define _TOPOLOGICAL_OPERATOR_STARBASE_DEFINE

/************************************************************************/
/* This file contains functions and data related to star cavity         */
/************************************************************************/
#include <VolumeMesh/Mesh/TetraMesh.h>
namespace VolumeMesh
{
	/*****************************************************************************/
	/*  orient()     Return a positive value if the point v lies below the plane */
	/*               passing through fp1, fp2, and fp3; "below" is defined so    */
	/*               that fp1, fp2, fp3 appear in counterclockwise order when    */
	/*               viewed from above the plane.  Returns a negative value if   */
	/*               v lies above the plane.  Returns zero if the points are     */
	/*               coplanar.                                                   */
	/*****************************************************************************/

	extern int orient(Point v, Point fp1, Point fp2, Point fp3);

	extern bool tetPoints(TetraMesh *mesh, TetraHandle handle, Point *p);
	extern bool facePoints(TetraMesh *mesh, HalfFaceHandle handle, Point *p);
	extern Point faceCenter(TetraMesh *mesh, TetraMesh::HalfFaceHandle fh);
	extern double tetVolume(Point p1, Point p2, Point p3, Point p4);
	extern void pointGrouping(int **neighbours, int *neighbourcnt, int pointcnt, int *pointcolour, int &colourcnt, float &time);
}
#endif