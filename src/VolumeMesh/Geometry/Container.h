/*---------------------------------------------------------------------------------------------------------------------
*  _               _  _             _   *
*   \             /  | \           / |  *     VolumeMesh : The Open Data Structure for Tetrahedral and Hexahedral Mesh
*    \           /   |  \         /  |  *     
*  	  \         /    |   \       /   |  *     Copyright(C) 2010 by Computer Aided Designed Group
*	   \       /     |    \     /    |  *     State Key Lab. of CAD & CG, Zhejiang University
*	    \     /      |     \   /     |  *
*	     \   /       |      \_/      |  *
*		  \_/        |               |_ * 
*                                       *
*-----------------------------------------------------------------------------------------------------------------------
* License
*
*    This file is part of VolumeMesh.
*
*    VolumeMesh is free software: you can redistribute it and/or modify       
*	it under the terms of the GNU Lesser General Public License as          
*	published by the Free Software Foundation, either version 3 of          
*	the License, or (at your option) any later version. 
*	
*	VolumeMesh distributed in the hope that it will be useful,            
*   but WITHOUT ANY WARRANTY; without even the implied warranty of      
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
*	GNU Lesser General Public License for more details.
*   This project is created by Chuhua Xian
*   Developers : Chuhua Xian,   chuhuaxian@gmail.com 
*                Xiaoshen Chen, chinimei@163.com
*
/---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*
*                                                                                                                      *
*                               VolumeMesh                                                                             *
*      Copyright (C) 2010 by Computer Aided Designed Group                                                             *
*        State Key Lab. of CAD & CG, Zhejiang University                                                               *
*         This project is created by Chuhua Xian, 2010.2                                                               *
*                     Email: chuhuaxian@gmail.com																	   *	
*         Modified by Xiaoshen Chen, 2010.03																		   *
*					  Email: chinimei@163.com																		   *
*----------------------------------------------------------------------------------------------------------------------*/ 

#ifndef _VOLUME_MESH_CONTAINER_H_
#define _VOLUME_MESH_CONTAINER_H_

//---------------------------------------------------------------------------------------------------------------------

#include <VolumeMesh/Mesh/Topology.h>
#include <VolumeMesh/Geometry/VMGeometry.h>
#include <VolumeMesh/Geometry/VMVector.h>
#include <vector>
#include <map>
#include <set>
#include <string>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	//typedef PointT<double> Point;
	typedef Vec3d Point;
	
	typedef std::vector<Tetrahedron>  TetraContainer;
	typedef std::vector<HalfFaceTH>   THFContainer;
	typedef std::vector<HalfEdgeTH>   THEContainer;
	typedef std::vector<Hexahedron>   HexadContainer;
	typedef std::vector<HalfFaceHH>   HHFContainer;
	typedef std::vector<HalfEdgeHH>   HHEContainer;
	typedef std::vector<Vertex>   VContainer;
	typedef std::vector<Point>    PContainer;	
	typedef std::vector<bool>     BRVContainer;     /**< boundary signs container  */
	typedef std::vector<Vec3d>    NContainer;       /**< normal container          */
	typedef std::set<PointHandle> BVContatiner;     /**< boundary vertex container */
	typedef std::map<std::string, int, std::less<std::string>> OPPContainer;

	typedef std::vector<VertexHandle> PVCContainer; /**< the container of the connections of PointHandle and VertexHandle */
	typedef std::map<PointHandle, std::vector<VertexHandle>> PVMContainer; /**< the container of the mapping from PointHandle to VertexHandle */

//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------

#endif