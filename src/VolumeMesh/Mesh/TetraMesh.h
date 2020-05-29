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
*					   Email: chinimei@163.com																		   *
*---------------------------------------------------------------------------------------------------------------------*/ 
/*---------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		       *
*					   Email: chuhuaxian@gmail.com																       *
*---------------------------------------------------------------------------------------------------------------------*/ 
/*---------------------------------------------------------------------------------------------------------------------
* Modified Information                                                                                                *
* Modified by Chuhua Xian                                                                                             *
* Email : chuhuaxian@gmail.com                                                                                        *
* Modified Date : 2010.8.17.                                                                                          *
*--------------------------------------------------------------------------------------------------------------------*/


#ifndef _VOLUME_MESH_TETRA_MESH_H_
#define _VOLUME_MESH_TETRA_MESH_H_
//---------------------------------------------------------------------------------------------------------------------

#include "VolumeMesh/Mesh/BaseMesh.h"
#include "VolumeMesh/Geometry/Container.h"
#include "VolumeMesh/Mesh/BaseMesh.h"
#include "VolumeMesh/Mesh/Iterators.h"
#include "VolumeMesh/Mesh/CirculatorT.h"
#include <algorithm>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	/// Tetra Main Control Class
	/** You can access, add, erase and iterate all the element of 
	*	the Tetra Topology.
	*/
	const int vhecTH_[4] = {3,0,1,2};  // container of first half edge relates to vertex
	//int vhfcTH_[4] = {1,2,3,0};  // container of first half face relates to vertex
	class TetraMesh : public BaseMesh
	{
	public:
		typedef THEContainer HEContainer;
		typedef THFContainer HFContainer;

		typedef Circulator::HedronHalfFaceIterT<TetraMesh>     HedronFaceIter;
		typedef Circulator::HalfFaceVertexIterT<TetraMesh>     HalfFaceVertexIter;
		typedef Circulator::HedronVertexIterT<TetraMesh>       HedronVertexIter;
		typedef Circulator::HedronHedronIterT<TetraMesh>       HedronHedronIter;
		typedef Circulator::VertexHedronIterT<TetraMesh>       VertexHedronIter;
		typedef Circulator::PointHedronIterT<TetraMesh>        PointHedronIter;    
		typedef Circulator::FaceHalfedgeIterT<TetraMesh>       FaceHalfedgeIter;
		typedef Circulator::HedronHalfEdgeIterT<TetraMesh>     HedronHalfEdgeIter;
		typedef Circulator::VertexHalfEdgeIterT<TetraMesh>     VertexHalfEdgeIter;
		typedef Circulator::VertexHalfFaceIterT<TetraMesh>     VertexHalfFaceIter;

	public:
		TetraMesh();
		~TetraMesh();
		TetraMesh(TetraMesh & _mesh) : BaseMesh(0x0001)
		{
			//this->_type = _mesh._type;
			this->_pc = _mesh._pc;
			this->_vc = _mesh._vc;
			this->_npoint = _mesh._npoint;
			this->_npolyhedron = _mesh._npolyhedron;
			this->_brv = _mesh._brv;
			this->_fn = _mesh._fn;
			this->_hec = _mesh._hec;
			this->_hfc = _mesh._hfc;
			this->_ophfc = _mesh._ophfc;
			this->_pvc = _mesh._pvc;
			this->_pvm = _mesh._pvm;
			this->_rfn = _mesh._rfn;
			this->_tetc = _mesh._tetc;
			
		}
	//-----------------------------------------------------------------------------------------------------------------
	public:
		// translate the handle to the item in the container
		Tetrahedron  & handle_to_entity(TetraHandle    & handle);
		HalfFaceTH   & handle_to_entity(HalfFaceHandle & handle);
		HalfEdgeTH   & handle_to_entity(HalfEdgeHandle & handle);
		Vertex       & handle_to_entity(VertexHandle   & handle);
		Point        & handle_to_entity(PointHandle    & handle);

		const Tetrahedron  & handle_to_entity(TetraHandle    & handle) const;
		const HalfFaceTH   & handle_to_entity(HalfFaceHandle & handle) const;
		const HalfEdgeTH   & handle_to_entity(HalfEdgeHandle & handle) const;
		const Vertex       & handle_to_entity(VertexHandle   & handle) const;
		const Point        & handle_to_entity(PointHandle    & handle) const;

		Tetrahedron & iter_to_entity(HedronIter & iter);
		//HalfFaceTH   & iter_to_entity(HalfFaceHandle & handle);
		//HalfEdgeTH   & iter_to_entity(HalfEdgeHandle & handle);
		Vertex      & iter_to_entity(VertexIter& iter);
		Point       & iter_to_entity(PointIter& iter);
	
		//Next Prev Opposite Radial operator
		Tetrahedron & next_tetrahedron(TetraHandle  & handle);
		HalfFaceTH  & next_half_face(HalfFaceHandle & handle);
		HalfEdgeTH  & next_half_edge(HalfEdgeHandle & handle);
		Tetrahedron & next_tetrahedron(Tetrahedron  & th_);
		HalfFaceTH  & next_half_face(HalfFaceTH     & hf_);
		HalfEdgeTH  & next_half_edge(HalfEdgeTH     & he_);

		TetraHandle    next_tetrahedron_handle(TetraHandle  & handle);
		HalfFaceHandle next_half_face_handle(HalfFaceHandle & handle);
		HalfEdgeHandle next_half_edge_handle(HalfEdgeHandle & handle);
		TetraHandle    next_tetrahedron_handle(Tetrahedron  & th_);
		HalfFaceHandle next_half_face_handle(HalfFaceTH     & hf_);
		HalfEdgeHandle next_half_edge_handle(HalfEdgeTH     & he_);

		const Tetrahedron & next_tetrahedron(TetraHandle  & handle) const;
		const HalfFaceTH  & next_half_face(HalfFaceHandle & handle) const;
		const HalfEdgeTH  & next_half_edge(HalfEdgeHandle & handle) const;
		const Tetrahedron & next_tetrahedron(Tetrahedron  & th_)    const;
		const HalfFaceTH  & next_half_face(HalfFaceTH     & hf_)    const;
		const HalfEdgeTH  & next_half_edge(HalfEdgeTH     & he_)    const;

		Tetrahedron & prev_tetrahedron(TetraHandle  & handle);
		HalfFaceTH  & prev_half_face(HalfFaceHandle & handle);
		HalfEdgeTH  & prev_half_edge(HalfEdgeHandle & handle);
		Tetrahedron & prev_tetrahedron(Tetrahedron  & th_);
		HalfFaceTH  & prev_half_face(HalfFaceTH     & hf_);
		HalfEdgeTH  & prev_half_edge(HalfEdgeTH     & he_);

		TetraHandle    prev_tetrahedron_handle(TetraHandle  & handle);
		HalfFaceHandle prev_half_face_handle(HalfFaceHandle & handle);
		HalfEdgeHandle prev_half_edge_handle(HalfEdgeHandle & handle);
		TetraHandle    prev_tetrahedron_handle(Tetrahedron  & th_);
		HalfFaceHandle prev_half_face_handle(HalfFaceTH     & hf_);
		HalfEdgeHandle prev_half_edge_handle(HalfEdgeTH     & he_);

		const Tetrahedron & prev_tetrahedron(TetraHandle  & handle) const; 
		const HalfFaceTH  & prev_half_face(HalfFaceHandle & handle) const;
		const HalfEdgeTH  & prev_half_edge(HalfEdgeHandle & handle) const;
		const Tetrahedron & prev_tetrahedron(Tetrahedron  & th_)    const;
		const HalfFaceTH  & prev_half_face(HalfFaceTH     & hf_)    const;
		const HalfEdgeTH  & prev_half_edge(HalfEdgeTH     & he_)    const;

		HalfEdgeTH  & radial_half_edge(HalfEdgeHandle   & handle);
		HalfEdgeTH  & mate_half_edge(HalfEdgeHandle     & handle);
		HalfFaceTH  & opposite_half_face(HalfFaceHandle & handle);
		Tetrahedron & get_tetrahedron(VertexHandle      & handle);

		HalfEdgeHandle radial_half_edge_handle(HalfEdgeHandle   & handle);
		HalfEdgeHandle mate_half_edge_handle(HalfEdgeHandle     & handle);
		HalfFaceHandle opposite_half_face_handle(HalfFaceHandle & handle);
		TetraHandle    tetrahedron_handle(VertexHandle          & handle);

		bool has_radial_half_edge(HalfEdgeHandle   &handle);
		bool has_opposite_half_face(HalfFaceHandle & handle);
		bool has_vertex(HalfFaceHandle hfHandle, VertexHandle vHandle);
		bool has_point(HalfFaceHandle hfHandle, PointHandle pHandle);
		bool has_boundary_face(PolyhedronHandle & hh_);
		bool has_boundary_face(HedronIter & hIter_);

		/// is boundary
		bool is_boundary(HalfFaceHandle   & hfh_);
		bool is_boundary(VertexHandle     & vh_);
		bool is_boundary(PolyhedronHandle & hh_);
		bool is_boundary(HedronHandle & hh_, unsigned int type_);
		bool is_boundary(HalfEdgeHandle &heh_);
		bool is_boundary(PointHandle &ph_);

		/// is reverse tetra
		bool is_reverse(TetraHandle & th);

		/// Short-cut from HalfEdgeTH to point
		PointHandle  start_point_handle(HalfEdgeHandle& handle);
		VertexHandle from_vertex_handle(HalfEdgeHandle hh_);
		VertexHandle to_vertex_handle(HalfEdgeHandle  hh_);


		/// clear 
		void release_mesh();

		/// add Tetra
		TetraHandle add_tetrahedron(Point &p1, Point &p2, Point &p3, Point &p4);

		/// add point
		PointHandle add_point(Point p);

		/// delete point
		int erase_point(PointHandle &ph);


		/// delete Tetra
		TetraHandle erase_tetrahedron(TetraHandle &handle);

		//int new_erase_tetrahedron(TetraHandle &handle);

		/// advanced operation
		int vertex_star(VertexHandle vh, std::vector<TetraHandle> &tetraVec);
		int point_star(PointHandle ph, std::vector<TetraHandle> &tetraVec);
		int edge_star(HalfEdgeHandle eh, std::vector<TetraHandle> &tetraVec);
		int vertex_half_edge_star(VertexHandle vh, std::vector<HalfEdgeHandle> &halfEdgeVec);
		int vertex_half_face_star(VertexHandle vh, std::vector<HalfFaceHandle> &halfFaceVec);

		HalfEdgeHandle & half_edge_handle(HalfFaceHandle & fh, VertexHandle & vh);
		HalfFaceHandle & half_face_handle(HalfEdgeHandle & handle);

		/// build the topology
		bool build_topology();
		bool build_opposite_half_face();

		/// topology operation : edge contract
		int edge_contract(HalfEdgeHandle &he_, PointHandle &ph_new);

		/// operation recover
		int recover_edge_contract();
	//-----------------------------------------------------------------------------------------------------------------
	public:
		Point & start_point(HalfEdgeHandle& he_);
		Point & point(PointHandle & handle);
		Point & point(VertexHandle & handle);
		Point & point(PointIter  & iter);
		Point & point(VertexIter & iter);

		HalfFaceHandle first_half_face_handle(TetraHandle th_);
		HalfFaceHandle first_half_face_handle(HedronIter tIter_);

		HalfEdgeHandle first_half_edge_handle(HalfFaceHandle fh_);

		HalfEdgeHandle first_hedron_half_edge_handle(HedronHandle & hh_);
		HalfEdgeHandle first_hedron_half_edge_handle(HedronIter   & hIter_);
		HalfFaceHandle first_vertex_half_face_handle(VertexHandle & vh_);
		HalfFaceHandle first_vertex_half_face_handle(VertexIter   & vIter_);
		HalfEdgeHandle first_vertex_half_edge_handle(VertexHandle & vh_);
		HalfEdgeHandle first_vertex_half_edge_handle(VertexIter   & vIter_);

		HalfEdgeHandle next_hedron_half_edge_handle(HedronHandle & hh_,  HalfEdgeHandle & heh_);
		HalfEdgeHandle next_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_);
		HalfEdgeHandle prev_hedron_half_edge_handle(HedronHandle & hh_,  HalfEdgeHandle & heh_);
		HalfEdgeHandle prev_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_);
		HalfFaceHandle next_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_);
		HalfFaceHandle next_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_);
		HalfFaceHandle prev_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_);
		HalfFaceHandle prev_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_);

		VertexHandle first_vertex_handle(HalfFaceHandle fh_);
		VertexHandle first_vertex_handle(PolyhedronHandle hh_);
		VertexHandle first_vertex_handle(HedronIter tIter_);
		VertexHandle next_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_);
		VertexHandle next_vertex_handle(HedronIter tIter_, VertexHandle cvh_);
		VertexHandle prev_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_);
		VertexHandle prev_vertex_handle(HedronIter tIter_, VertexHandle cvh_);

		void update_face_normals();
		Normal update_face_normal(HalfFaceHandle & hfh_);

		Normal normal(HalfFaceHandle & hfh_);
		Normal normal(const HalfFaceHandle & hfh_) const;

		void update_boundary_point();
	//-----------------------------------------------------------------------------------------------------------------
	public:
		/** adjust the order of the vertices of the tetrahedron, make the normal of the face pointing outside
		*   \param in hIter_ the iterator of the tetrahedron
		*/
		// need testing
		void adjust_tetrahedron(HedronIter & hIter_)
		{
			unsigned int i;
			i = hIter_.handle().idx();
			Normal n[2];
			n[0] = (_pc[_vc[i * 4 + 2].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]) % 
				   (_pc[_vc[i * 4 + 3].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]);
			n[1] = _pc[_vc[i * 4 + 1].point_handle().idx()] - _pc[_vc[i * 4].point_handle().idx()];

			VertexHandle curr;
			if ((n[0] | n[1]) < 0)
			{
				// adjust _pvc loop
				curr = _pvc[_vc[i*4+1].point_handle().idx()];
				if (curr != _vc[i*4+1].handle())
				{
					while (handle_to_entity(curr).next_vertex_handle() != _vc[i*4+1].handle())
						curr = handle_to_entity(curr).next_vertex_handle();
					assert(curr!=VertexHandle(-1));
					handle_to_entity(curr).set_next_vertex_handle(handle_to_entity(handle_to_entity(curr).next_vertex_handle()).next_vertex_handle());
				}
				else
					_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();

				curr = _pvc[_vc[i*4+3].point_handle().idx()];
				if (curr != _vc[i*4+3].handle())
				{
					while (handle_to_entity(curr).next_vertex_handle() != _vc[i*4+3].handle())
						curr = handle_to_entity(curr).next_vertex_handle();
					assert(curr!=VertexHandle(-1));
					handle_to_entity(curr).set_next_vertex_handle(handle_to_entity(handle_to_entity(curr).next_vertex_handle()).next_vertex_handle());
				}
				else
					_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();

				// swap & reset two vertex
				PointHandle tph = _vc[i*4+1].point_handle();
				_vc[i*4+1].set_point_handle(_vc[i*4+3].point_handle());
				_vc[i*4+3].set_point_handle(tph);
				_vc[i*4+1].set_next_vertex_handle(_pvc[_vc[i*4+1].point_handle().idx()]);
				_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].handle();
				_vc[i*4+3].set_next_vertex_handle(_pvc[_vc[i*4+3].point_handle().idx()]);
				_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].handle();
			}		
		}

		/** adjust the order of the vertices of the tetrahedron, make the normal of the faces pointing outside
		*   \param in hh_ the handle of the tetrahedron
		*/
		void adjust_tetrahedron(HedronHandle & hh_)
		{
			unsigned int i;
			i = hh_.idx();
			Normal n[2];
			n[0] = (_pc[_vc[i * 4 + 2].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]) % 
				   (_pc[_vc[i * 4 + 3].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]);
			n[1] = _pc[_vc[i * 4 + 1].point_handle().idx()] - _pc[_vc[i * 4].point_handle().idx()];

			VertexHandle curr;
			if ((n[0] | n[1]) < 0)
			{
				// adjust _pvc loop
				curr = _pvc[_vc[i*4+1].point_handle().idx()];
				if (curr != _vc[i*4+1].handle())
				{
					while (handle_to_entity(curr).next_vertex_handle() != _vc[i*4+1].handle())
						curr = handle_to_entity(curr).next_vertex_handle();
					assert(curr!=VertexHandle(-1));
					handle_to_entity(curr).set_next_vertex_handle(handle_to_entity(handle_to_entity(curr).next_vertex_handle()).next_vertex_handle());
				}
				else
					_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();

				curr = _pvc[_vc[i*4+3].point_handle().idx()];
				if (curr != _vc[i*4+3].handle())
				{
					while (handle_to_entity(curr).next_vertex_handle() != _vc[i*4+3].handle())
						curr = handle_to_entity(curr).next_vertex_handle();
					assert(curr!=VertexHandle(-1));
					handle_to_entity(curr).set_next_vertex_handle(handle_to_entity(handle_to_entity(curr).next_vertex_handle()).next_vertex_handle());
				}
				else
					_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();

				// swap & reset two vertex
				PointHandle tph = _vc[i*4+1].point_handle();
				_vc[i*4+1].set_point_handle(_vc[i*4+3].point_handle());
				_vc[i*4+3].set_point_handle(tph);
				_vc[i*4+1].set_next_vertex_handle(_pvc[_vc[i*4+1].point_handle().idx()]);
				_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].handle();
				_vc[i*4+3].set_next_vertex_handle(_pvc[_vc[i*4+3].point_handle().idx()]);
				_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].handle();
			}		
		}


	    /** set the geometry position of a point handle
		*   \param in _ph the vertex handle
		*   \param in _p the position of the point
		*/
		void set_point(PointHandle _ph, Point _p)
		{
			_pc[_ph.idx()] = _p;
		}

		/** set the point of a vertex handle
		*/
		void set_point(VertexHandle _vh, Point & _p)
		{
			_pc[handle_to_entity(_vh).point_handle().idx()] = _p;
		}

		/** set the point of a vertex handle
		*/
		void set_point(VertexIter _vIter, Point & _p)
		{
			_pc[handle_to_entity(_vIter.handle()).point_handle().idx()] = _p;
		}

		/** retrieve the point handle from the vertex handle
		*/
		PointHandle point_handle(VertexHandle _vh)
		{
			return handle_to_entity(_vh).point_handle();
		}
		
		/** retrieve the point handle from the vertex handle
		*/
		VertexHandle vertex_handle(PointHandle _vh)
		{
			return VertexHandle(_pvc[_vh.idx()]);
		}

		/** calculate the mate dihedral angle of two mate half face
		*   \param in heh_ the half edge of one face
		*   \return the dihedral of the two half face
		*/
		Scalar calc_mate_dihedral_angle(HalfEdgeHandle heh_)
		{
			Point v[4];
			v[0] = point(from_vertex_handle(heh_));
			v[1] = point(to_vertex_handle(heh_));
			v[2] = point(to_vertex_handle(next_half_edge_handle(heh_)));

			v[3] = point(to_vertex_handle(next_half_edge_handle(mate_half_edge_handle(heh_))));

			Normal n[2];

			n[0] = (v[1] - v[0]) % (v[2] - v[0]);
			n[0].normalize();
			n[1] = (v[3] - v[0]) % (v[1] - v[0]);
			n[1].normalize();

			Scalar da_cos;
			da_cos = acos((std::max)(-1.0, (std::min)(1.0, n[0] | n[1])));
			Scalar da_sign;
			da_sign = (n[0] % n[1]) | (v[1] - v[0]);
			if (da_sign > 0)
			{
				return M_PI + da_cos;
			} 
			return M_PI - da_cos;
		}

		/** calculate the radial dihedral angle of two mate half face
		*   \param in heh_ the half edge of one face
		*   \return the dihedral of the two half face
		*/
		Scalar calc_radial_dihedral_angle(HalfEdgeHandle heh_)
		{
			Point v[4];
			v[0] = point(from_vertex_handle(heh_));
			v[1] = point(to_vertex_handle(heh_));
			v[2] = point(to_vertex_handle(next_half_edge_handle(heh_)));

			v[3] = point(to_vertex_handle(next_half_edge_handle(radial_half_edge_handle(heh_))));

			Normal n[2];

			n[0] = (v[1] - v[0]) % (v[2] - v[0]);
			n[0].normalize();
			n[1] = (v[3] - v[0]) % (v[1] - v[0]);
			n[1].normalize();

			Scalar da_cos;
			da_cos = acos((std::max)(-1.0, (std::min)(1.0, n[0] | n[1])));
			Scalar da_sign;
			da_sign = (n[0] % n[1]) | (v[1] - v[0]);
			if (da_sign > 0)
			{
				return M_PI + da_cos;
			} 
			return da_cos;
		}

		/** retrieve the boundary hedron handle though the point handle of a boundary half face determined by three points
		*   if the half face determined by the three points is not a boundary half face, then return a invalid handle
		*/
		HedronHandle boundary_half_face_hedron_handle(PointHandle & ph0_, PointHandle & ph1_, PointHandle & ph2_)
		{
			char buf[128];
			int pointHandle[3];			
			pointHandle[0] = ph0_.idx();
			pointHandle[1] = ph1_.idx();
			pointHandle[2] = ph2_.idx();
			std::sort(pointHandle, pointHandle + 3);
			sprintf_s(buf, "%d %d %d", pointHandle[0], pointHandle[1], pointHandle[2]);
			if (_ophfc.find(buf) != _ophfc.end())
			{
				HalfFaceHandle hfh(_ophfc[buf]);
				if ((hfh.idx() < 0) || (hfh.idx() >= (int)_hfc.size()))
				{
					return HedronHandle(-1);
				}
				return handle_to_entity(hfh).hedron_handle();
			}
			return HedronHandle(-1);
		}

	//-----------------------------------------------------------------------------------------------------------------

	public:
		/** return the first iterator of the HalfFaceIter of the tetrahedron
		 *  \param th_ the handle of the tetrahedron
		 *  \return HedronFaceIter the iterator.
		 */
		HedronFaceIter hedron_face_iter(TetraHandle th_)
		{
			return HedronFaceIter(*this, first_half_face_handle(th_));
		}
		/** return the first iterator of the HalfFaceIter of the tetrahedron
		*  \param hIter_ the iterator of the tetrahedron
		*  \return HedronFaceIter the iterator.
		*/
		HedronFaceIter hedron_face_iter(HedronIter hIter_)
		{
			return HedronFaceIter(*this, first_half_face_handle(hIter_.handle()));
		}


		/** return the first iterator of the HalfFaceVertexIter of the half face
		 *  \param hfh_ the handle of the half face
		 *  \return HalfFaceVertexIter the iterator
		 */
		HalfFaceVertexIter half_face_vertex_iter(HalfFaceHandle hfh_)
		{
			return HalfFaceVertexIter(*this, first_half_edge_handle(hfh_));
		}

		/** return the first iterator of the HedronVertexIter of the tetrahedron
		*  \param th_ the handle of the hedron
		*  \return HedronVertexIter the iterator of the vertex
		*/
		HedronVertexIter hedron_vertex_iter(PolyhedronHandle th_)
		{
			return HedronVertexIter(*this, th_, first_vertex_handle(th_));
		}

		/** return the first iterator of the HedronVertexIter of the tetrahedron
		*  \param tIter_ the handle of the hedron
		*  \return HedronVertexIter the iterator of the vertex
		*/
		HedronVertexIter hedron_vertex_iter(HedronIter tIter_)
		{
			return HedronVertexIter(*this, tIter_.handle(), first_vertex_handle(tIter_));
		}

		/** return the first iterator of HedronHedronIter the polyhedron 
		*   \param ph_ the handle of the polyhedron
		*   \return HedronHedronIter the iterator of the polyhedron
		*/
		HedronHedronIter hedron_hedron_iter(PolyhedronHandle ph_)
		{
			return HedronHedronIter(*this, first_half_face_handle(ph_));
		}

		/** return the first iterator of VertexHedronIter the polyhedron 
		*   \param vh_ the handle of the vertex
		*   \return VertexHedronIter 
		*/
		VertexHedronIter vertex_hedron_iter(const VertexHandle & vh_)
		{
			std::vector<HedronHandle> hedrons;
			vertex_star(vh_, hedrons);
			return VertexHedronIter(*this, hedrons);
		}
		
		/** return the first iterator of VertexHalfEdgeIter the half edge 
		*   \param vh_ the handle of the vertex
		*   \return VertexHalfEdgeIter 
		*/
		VertexHalfEdgeIter vertex_half_edge_iter(VertexHandle & vh_)
		{
			return VertexHalfEdgeIter(*this, vh_, first_vertex_half_edge_handle(vh_));
		}

		/** return the first iterator of VertexHalfFaceIter the half face
		*   \param vh_ the handle of the vertex
		*   \return VertexHalfFaceIter 
		*/
		VertexHalfFaceIter vertex_half_face_iter(VertexHandle & vh_)
		{
			return VertexHalfFaceIter(*this, vh_, first_vertex_half_face_handle(vh_));
		}

		/** return the first iterator of the PointHedronIter of the point handle
		*   \param ph_ the handle of the point
		*   \return the first iterator of the PointHedronIter
		*/
		PointHedronIter point_hedron_iter(const PointHandle & ph_)
		{
			std::vector<HedronHandle> hedrons;
			vertex_star(_pvc[ph_.idx()], hedrons);
			return PointHedronIter(*this, hedrons);
		}

		/** return the first iterator of VertexHedronIter the polyhedron 
		*   \param vh_ the handle of the vertex
		*   \return VertexHedronIter 
		*/
		VertexHedronIter vertex_hedron_iter(VertexIter vIter_)
		{
			std::vector<HedronHandle> hedrons;
			vertex_star(vIter_.handle(), hedrons);
			return VertexHedronIter(*this, hedrons);
		}

		/** return the first iterator of the PointHedronIter of the point handle
		*   \param ph_ the handle of the point
		*   \return the first iterator of the PointHedronIter
		*/
		PointHedronIter point_hedron_iter(PointIter pIter_)
		{
			std::vector<HedronHandle> hedrons;
			vertex_star(_pvc[pIter_.handle().idx()], hedrons);
			return PointHedronIter(*this, hedrons);
		}

		/** return the first FaceHalfedgeIter of the half face handle
		*   \param hf_ the handle of the half face
		*   \return the first iterator of the FaceHalfedgeIter
		*/
		FaceHalfedgeIter face_half_edge_iter(HalfFaceHandle hf_)
		{
			return FaceHalfedgeIter(*this, first_half_edge_handle(hf_));
		}

		/** return the first FaceHalfedgeIter of the half face handle
		*   \param hf_ the handle of the half face
		*   \return the first iterator of the FaceHalfedgeIter
		*/
		FaceHalfedgeIter fhe_iter(HalfFaceHandle hf_)
		{
			return FaceHalfedgeIter(*this, first_half_edge_handle(hf_));
		}

		/** return the first half edge of a hedron
		*   \param in hh_ the hedron iter
		*   \return the first half edge iter
		*/
		HedronHalfEdgeIter hedron_half_edge_iter(HedronIter & hIter_)
		{
			return HedronHalfEdgeIter(*this, hIter_.handle());
		}

		/** return the first half edge handle of a hedron
		*   \param in hh_ the hedron handle
		*   \return the first half edge iter
		*/
		HedronHalfEdgeIter hedron_half_edge_iter(HedronHandle & hh_)
		{
			return HedronHalfEdgeIter(*this, hh_);
		}

		
	//-----------------------------------------------------------------------------------------------------------------
	public:
		HedronIter hedrons_begin()
		{
			return HedronIter(first_valid_tetrahedron_handle().idx());
		}

		HedronIter hedrons_end()
		{
			return HedronIter(n_polyhedron());
		}

		/** return the hedron handle from the half face
		*   \param hfh the half face handle
		*   \return PolyhedronHandle the handle of the Polyhedron
		*/
		PolyhedronHandle hedron_handle(HalfFaceHandle hfh)
		{
			if ((hfh.idx() < 0) || (hfh.idx() >= (int)_hfc.size()))
			{
				return PolyhedronHandle(-1);
			}
			return handle_to_entity(hfh).hedron_handle();
		}

		/** return the first adjacent hedron handle from the current hedron handle
		*   \param ph the polyhedron handle
		*   \return PolyhedronHandle the handle of the Polyhedron
		*/
		PolyhedronHandle first_adjacent_hedron_handle(PolyhedronHandle ph)
		{
			if ((ph.idx() < 0) || (ph.idx() >= (int)_tetc.size()))
			{
				return PolyhedronHandle(-1);
			}
			return handle_to_entity(handle_to_entity(ph).first_half_face_handle()).hedron_handle();
		}	

		/** retrieve the adjacent half face handle of two hedrons
		*   \return the half face handle of the first hedron
		*/
		HalfFaceHandle adjacent_half_face_handle(HedronHandle h0_, HedronHandle h1_)
		{
			std::set<HalfFaceHandle> hfs;
			for (HedronFaceIter hf_it = hedron_face_iter(h1_); hf_it ; ++ hf_it)
			{
				hfs.insert(opposite_half_face_handle(hf_it.handle()));
			}

			for (HedronFaceIter hf_it = hedron_face_iter(h0_); hf_it ; ++ hf_it)
			{
				if (hfs.find(hf_it.handle()) != hfs.end())
				{
					return hf_it.handle();
				}
			}
			return HalfFaceHandle(-1);
		}

		int size_tetrahedron()
		{
			return _tetc.size();
		}

		int size_tetrahedron() const
		{
			return _tetc.size();
		}

		int size_halfface()
		{
			return _hfc.size();
		}

		int size_halfedge()
		{
			return _hec.size();
		}

		HalfFaceHandle vertex_opposite_half_face_handle(VertexHandle vh)
		{
			HalfFaceTH hff = _hfc[vh];
			return hff.handle();
		}

		HalfEdgeHandle vertex_first_half_edge_handle(VertexHandle vh)
		{
			return HalfEdgeHandle((vh>>2)*12+vhecTH_[vh%4]);
		}
		//-----------------------------------------------------------------------------------------------------------------
		//API about _free_list

		bool is_valid(TetraHandle &th_)
		{
			// valid handle and not in the free list
			if ( th_>=0 && th_<size_tetrahedron() && find(_free_tetra_list.begin(), _free_tetra_list.end(), th_)==_free_tetra_list.end() )
				return true;
			else
				return false;
		}

		bool is_valid(TetraHandle &th_) const
		{
			// valid handle and not in the free list
			if ( th_>=0 && th_<size_tetrahedron() && find(_free_tetra_list.begin(), _free_tetra_list.end(), th_)==_free_tetra_list.end() )
				return true;
			else
				return false;
		}

		bool is_valid(HalfFaceHandle &hfh_)
		{
			return is_valid(this->handle_to_entity(hfh_).hedron_handle());
		}

		bool is_valid(HalfEdgeHandle &heh_)
		{
			return is_valid(this->handle_to_entity(heh_).hedron_handle());
		}

		bool is_valid(VertexHandle &vh_)
		{
			return is_valid(this->tetrahedron_handle(vh_));
		}

		bool is_valid(PointHandle &ph_)
		{
			// valid handle and not in the free list
			if ( ph_>=0 && ph_<size_point() && find(_free_point_list.begin(), _free_point_list.end(), ph_)==_free_point_list.end() )
				return true;
			else
				return false;
		}

		//return next valid handle	
		//return -1 if no handle at all
		TetraHandle next_valid_tetrahedon_handle(TetraHandle &th_)
		{
			//assert if input TetraHandle is out of range;
			assert ( (th_.idx() < _npolyhedron ) && (th_.idx() >= 0) );
			int idx = th_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==_npolyhedron-1 )
					idx = 0;
				else 
					idx++;
				//find out valid tetrahedeon
				if( is_valid(TetraHandle(idx)) )
					return TetraHandle(idx);
			}
			while(idx!=th_.idx());
			
			return TetraHandle(-1);
		}
		
		TetraHandle next_valid_tetrahedon_handle(TetraHandle &th_) const
		{
			//assert if input TetraHandle is out of range;
			assert ( (th_.idx() < _npolyhedron ) && (th_.idx() >= 0) );
			int idx = th_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==_npolyhedron-1 )
					idx = 0;
				else 
					idx++;
				//find out valid tetrahedeon
				if( is_valid(TetraHandle(idx)) )
					return TetraHandle(idx);
			}
			while(idx!=th_.idx());

			return TetraHandle(-1);
		}

		//return prev valid handle
		//return -1 if no handle at all
		TetraHandle prev_valid_tetrahedon_handle(TetraHandle &th_)
		{
			//assert if input TetraHandle is out of range;
			assert ( (th_.idx() < _npolyhedron ) && (th_.idx() >= 0) );
			int idx = th_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==0 )
					idx = _npolyhedron-1;
				else 
					idx--;
				//find out valid tetrahedeon
				if( is_valid(TetraHandle(idx)) )
					return TetraHandle(idx);
			}
			while(idx!=th_.idx());

			return TetraHandle(-1);
		}

		TetraHandle prev_valid_tetrahedon_handle(TetraHandle &th_) const
		{
			//assert if input TetraHandle is out of range;
			assert ( (th_.idx() < _npolyhedron ) && (th_.idx() >= 0) );
			int idx = th_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==0 )
					idx = _npolyhedron-1;
				else 
					idx--;
				//find out valid tetrahedeon
				if( is_valid(TetraHandle(idx)) )
					return TetraHandle(idx);
			}
			while(idx!=th_.idx());

			return TetraHandle(-1);
		}

		//return first valid tetrahedron handle
		TetraHandle first_valid_tetrahedron_handle()
		{
			return next_valid_tetrahedon_handle(TetraHandle(_npolyhedron-1));
		}

		TetraHandle last_valid_tetrahedron_handle()
		{
			return prev_valid_tetrahedon_handle(TetraHandle(0));
		}

		void add_to_free_tetra_list(TetraHandle &th_)
		{
			if (find(_free_tetra_list.begin(), _free_tetra_list.end(), th_)==_free_tetra_list.end())
				_free_tetra_list.push_back(th_);
		}

		void add_to_free_point_list(PointHandle &ph_)
		{
			if (find(_free_point_list.begin(), _free_point_list.end(), ph_)==_free_point_list.end())
				_free_point_list.push_back(ph_);
		}

		//Get first free tetra location, or return -1 if empty
		TetraHandle pop_from_free_tetra_list()
		{
			if( _free_tetra_list.size()>0 )
				return *_free_tetra_list.begin();
			else
				return TetraHandle(-1);
		}

		//Get first free tetra location, then erase it. or return -1 if empty
		TetraHandle pop_erase_from_free_tetra_list()
		{
			if ( _free_tetra_list.size()>0 )
			{
				TetraHandle th = *_free_tetra_list.begin();
				_free_tetra_list.erase(_free_tetra_list.begin());
				return th;
			}
			else
				return TetraHandle(-1);
		}

		PointHandle pop_from_free_point_list()
		{
			if( _free_point_list.size()>0 )
				return *_free_point_list.begin();
			else
				return PointHandle(-1);
		}

		PointHandle pop_erase_from_free_point_list()
		{
			if ( _free_point_list.size()>0 )
			{
				PointHandle ph = *_free_point_list.begin();
				_free_point_list.erase(_free_point_list.begin());
				return ph;
			}
			else
				return PointHandle(-1);
		}

		// clean _free_tetra_list and _free_point_list & filling the vacant
		void clean_garbage();
		// find whether there is a valid point in point vector 
		PointHandle find_point(PointFindIF<Point>& pif);

	//-----------------------------------------------------------------------------------------------------------------
	/**
	* members of the TetraMesh
	*/	
	protected:
		TetraContainer _tetc;
		HFContainer _hfc;
		HEContainer _hec;
		std::vector<TetraHandle> _free_tetra_list;
		std::vector<PointHandle> _free_point_list;

		/// operation recover data

		//-----------------------------------edge contract recover data--------------------------------//
		Point rec_ec_from_point;   /**< contracted half_edge's from_point  */
		Point rec_ec_to_point;     /**< contracted half_edge's to_point    */
		PointHandle rec_ec_from_point_handle;  /**< contracted half_edge's from_point_handle  */
		PointHandle rec_ec_to_point_handle;    /**< contracted half_edge's t_point_handle     */
		VertexHandle rec_ec_to_point_vertex_handle;     /**< vertexHandle in _pvc of to_point */
		std::vector<PointHandle> rec_ec_erase_point;    /**< erased pointHandles in           */
		std::vector<TetraHandle> rec_ec_erase_tetra;    /**< set of  erased tetrahedrons in edge contracting    */
		std::vector<VertexHandle> rec_ec_to_point_vertex_container;    /**< set of VertexHandles that their original PointHanle is to_point_handle     */
	};
	

//---------------------------------------------------------------------------------------------------------------------
} // namespace VolumeMesh

//---------------------------------------------------------------------------------------------------------------------
#endif // _VOLUME_MESH_TETRA_MESH_H_ defined
