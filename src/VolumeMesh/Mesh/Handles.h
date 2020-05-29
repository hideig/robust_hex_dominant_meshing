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
 *                                                                                                                     *
 *                               VolumeMesh                                                                            *
 *      Copyright (C) 2010 by Computer Aided Designed Group                                                            *
 *        State Key Lab. of CAD & CG, Zhejiang University                                                              *
 *         This project is created by Chuhua Xian, 2010.2                                                              *
 *                     Email: chuhuaxian@gmail.com																	   *	
 *         Modified by Xiaoshen Chen, 2010.03																		   *
 *					   Email: chinimei@163.com																		   *
 *---------------------------------------------------------------------------------------------------------------------*/ 
/*---------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		       *
*					   Email: chuhuaxian@gmail.com																       *
*---------------------------------------------------------------------------------------------------------------------*/ 



#ifndef _VOLUME_MESH_HANDLES_H_
#define _VOLUME_MESH_HANDLES_H_

#include <VolumeMesh/Mesh/BaseHandle.h>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	/// Handle for a vertex entity
	struct VertexHandle : public BaseHandle
	{
		explicit VertexHandle(int _idx = -1) : BaseHandle(_idx) 
		{
		}
	};

	/// Handle for a HalfEdgeTH entity
	struct HalfEdgeHandle : public BaseHandle
	{
		explicit HalfEdgeHandle(int _idx = -1) : BaseHandle(_idx) 
		{
		}
	};

	/// Handle for a edge entity
	struct EdgeHandle : public BaseHandle
	{
		explicit EdgeHandle(int _idx = -1) : BaseHandle(_idx) 
		{
		}
	};
	
	/// Handle for a HalfFace entity
	struct HalfFaceHandle : public BaseHandle
	{
		explicit HalfFaceHandle(int _idx = -1) : BaseHandle(_idx)
		{
		}
	};

	/// Handle for a face entity
	struct FaceHandle : public BaseHandle
	{
		explicit FaceHandle(int _idx = -1) : BaseHandle(_idx) 
		{
		}
	};

	/// Handle for a point
	struct PointHandle : public BaseHandle
	{
		explicit PointHandle(int _idx = -1) : BaseHandle(_idx)
		{
		}
	};

	/// Handle for the polyhedron 
	struct PolyhedronHandle : public BaseHandle
	{
		explicit PolyhedronHandle(int _idx = -1) : BaseHandle(_idx)
		{
		}
	};

	/// Handle for a volume entity
	//struct TetraHandle : public PolyhedronHandle
	//{
	//	explicit TetraHandle(int _idx = -1) : PolyhedronHandle(_idx)
	//	{
	//	}
	//	bool operator = (PolyhedronHandle & h_)
	//	{
	//		setIdx(h_.idx());
	//		return true;
	//	}
	//	bool operator == (PolyhedronHandle & h_) const
	//	{
	//		return h_.idx() == idx();
	//	}
	//};

	/// Handle for a volume entity
	//struct HexadHandle : public PolyhedronHandle
	//{
	//	explicit HexadHandle(int _idx = -1) : PolyhedronHandle(_idx)
	//	{
	//	}

	//	bool operator = (PolyhedronHandle & h_)
	//	{
	//		setIdx(h_.idx());
	//		return true;
	//	}

	//	bool operator == (PolyhedronHandle & h_) const
	//	{
	//		return h_.idx() == idx();
	//	}
	//};

//---------------------------------------------------------------------------------------------------------------------
} // namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------

#endif // _VOLUME_MESH_HANDLES_H_ defined
