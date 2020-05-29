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
 *                     Email: chuhuaxian@gmail.com                                                                     *
 *---------------------------------------------------------------------------------------------------------------------*/ 

/**********************************************************************************************************************
 * this file defines the BaseMesh class, it is the base class for all mesh                                            *
 * $ written by xianc $                                                                                               *
 * $ Version 0.1 $                                                                                                    *
 * $ Date: 2010-02-08 $                                                                                               *
 **********************************************************************************************************************/
/*---------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		       *
*					   Email: chuhuaxian@gmail.com																       *
*---------------------------------------------------------------------------------------------------------------------*/ 


#ifndef _VOLUME_MESH_BASE_MESH_H_
#define _VOLUME_MESH_BASE_MESH_H_
//---------------------------------------------------------------------------------------------------------------------

#include <VolumeMesh/System/config.h>
#include <VolumeMesh/Mesh/Handles.h>
#include <VolumeMesh/Mesh/Iterators.h>
#include <VolumeMesh/Geometry/VMVector.h>
#include <volumeMesh/Mesh/Topology.h>
#include <VolumeMesh/Geometry/Container.h>

#include <vector>

//---------------------------------------------------------------------------------------------------------------------


namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------
	class BaseMesh 
	{
	//-----------------------------------------------------------------------------------------------------------------
	//-------------the typedef area------------------------------------//
	public:
		typedef Vec3d            Normal;
		typedef Point            Point;
		typedef Point::Scalar    Scalar;
		typedef PolyhedronHandle HedronHandle;
		typedef VertexHandle     VertexHandle;
		typedef PointHandle      PointHandle;
		typedef HalfFaceHandle   HalfFaceHandle;
		typedef HalfEdgeHandle   HalfEdgeHandle;

		typedef VertexIterator     VertexIter;
		typedef PointIterator      PointIter;
		typedef PolyhedronIterator HedronIter;
	//-----------------------------------------------------------------------------------------------------------------
	public:
		BaseMesh();
		BaseMesh(unsigned int type_);
		virtual ~ BaseMesh();
	//-----------------------------------------------------------------------------------------------------------------
	public:
		/** Retrieve the type of the mesh 
		 *  \return unsigned int . 0x0001 for tetrahedral mesh, 0x0002 for hexahedral mesh
		 */
		unsigned int mesh_type()
		{
			return _type;
		}
	//-----------------------------------------------------------------------------------------------------------------

	//-------------------------------virtual functions---------------------------------------------------------------//
	public:
		virtual bool build_topology() = 0;
		virtual bool is_boundary(HedronHandle & hh_, unsigned int type_) = 0;
			
	//-----------------------------------------------------------------------------------------------------------------
	public:
		/** Retrieve the vertex container
		*   \return the vertex container
		*/
		const VContainer & vertex_container() const
		{
			return _vc;
		}

		/** set the vertex container
		*/
		void set_vertex_container(VContainer & vc_)
		{
			_vc = vc_;
		}

		/** Retrieve the point container
		*   \return the point container
		*/
		const PContainer & point_container() const
		{
			return _pc;
		}

		/** set the point container
		*/
		void set_point_container(PContainer & pc_)
		{
			_pc = pc_;
		}


		/** Retrieve the number of the points
		 * \return the number of the points
		 */
		int n_point()
		{
			return _npoint;
		}

		/** Retrieve the number of the polyhedrons
		 * \return the number of the polyhedrons
		 */
		int n_polyhedron()
		{
			return _npolyhedron;
		}

		/** request face normals
		*/
		void request_face_normals()
		{
			_rfn = true;
		}

		///** set the geometry position of a point handle
		//*   \param in _ph the vertex handle
		//*   \param in _p the position of the point
		//*/
		//void set_point(const PointHandle & _ph, const Point & _p)
		//{
		//	_pc[_ph.idx()] = _p;
		//}

		//void set_point(PointHandle _ph, Point _p)
		//{
		//	_pc[_ph.idx()] = _p;
		//}


	//-----------------------------------------------------------------------------------------------------------------

	public:

		VertexIter vertices_begin()
		{
			return VertexIter(0);
		}

		VertexIter vertices_end()
		{
			return VertexIter(_vc.size());
		}
		PointIter points_begin()
		{
			return PointIter(0);
		}

		PointIter points_end()
		{
			return PointIter(_pc.size());
		}
	//-----------------------------------------------------------------------------------------------------------------
	
		int size_point()
		{
			return _pc.size();
		}

	//-----------------------------------------------------------------------------------------------------------------
	protected:
		VContainer   _vc;     /**< the vertices container           */
		PContainer   _pc;	  /**< the points container             */
		OPPContainer _ophfc;  /**< the opposite half face container */

		PVCContainer _pvc;    /**< the container of the PointHandle and VertexHandle */
		BRVContainer _brv;    /**< the boundary signs for points    */
		PVMContainer _pvm;    /**< the container of the mapping from PointHandle to VertexHandle */
		                      /**< better check if pointHandle and vertexHandle is valid before using them */

	protected:
		NContainer _fn;       /**< the face normal container        */
		bool       _rfn;   

	protected:
		int _npolyhedron; /**< the number of the polyhedrons */
		int _npoint;      /**< the number of the points      */
	private:
		unsigned int _type;   /**< the type of the mesh, 0x0001 for tetrahedral mesh, 0x0002 for hexahedral mesh */
	};
	
//---------------------------------------------------------------------------------------------------------------------
} // namespace VolumeMesh

//---------------------------------------------------------------------------------------------------------------------
#endif