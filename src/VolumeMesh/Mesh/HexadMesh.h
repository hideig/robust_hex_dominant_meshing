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

#ifndef _VOLUME_MESH_HEXADMESH_H_
#define _VOLUME_MESH_HEXADMESH_H_

//---------------------------------------------------------------------------------------------------------------------

#include <VolumeMesh/Mesh/BaseMesh.h>
#include <VolumeMesh/Geometry/Container.h>
#include <VolumeMesh/Mesh/BaseMesh.h>
#include <VolumeMesh/Mesh/Iterators.h>
#include <VolumeMesh/Mesh/CirculatorT.h>
#include <algorithm>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	/// Hexad Main Control Class
	/** You can access, add, erase and iterate all the element of 
	*	the Hexad Topology.
	*/
	const int vhecHH_[8] = {12,8,0,1,4,5,3,2};  // container of first half edge relates to vertex
	//int vhfcHH_[24] = {3,4,5, 2,3,5, 0,2,3, 0,3,4, 1,4,5, 1,2,5, 0,1,2, 0,1,4};  // container of first half face relates to vertex
	class HexadMesh : public BaseMesh
	{
	public:
		typedef HHEContainer HEContainer;
		typedef HHFContainer HFContainer;		

		typedef Circulator::HedronHalfFaceIterT<HexadMesh>       HedronFaceIter;
		typedef Circulator::HalfFaceVertexIterT<HexadMesh>       HalfFaceVertexIter;
		typedef Circulator::HedronVertexIterT<HexadMesh>         HedronVertexIter;
		typedef Circulator::HedronHedronIterT<HexadMesh>         HedronHedronIter;
		typedef Circulator::VertexHedronIterT<HexadMesh>         VertexHedronIter;
		typedef Circulator::PointHedronIterT<HexadMesh>          PointHedronIter;
		typedef Circulator::FaceHalfedgeIterT<HexadMesh>         FaceHalfedgeIter;
		typedef Circulator::HedronHalfEdgeIterT<HexadMesh>       HedronHalfEdgeIter;
		typedef Circulator::VertexHalfEdgeIterT<HexadMesh>       VertexHalfEdgeIter;
		typedef Circulator::VertexHalfFaceIterT<HexadMesh>       VertexHalfFaceIter;

	public:
		HexadMesh() : BaseMesh(0x0002)
		{
		}

		HexadMesh(HexadMesh & _mesh) : BaseMesh(0x0002)
		{
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
			this->_rfn = _mesh._rfn;
			this->_hexac = _mesh._hexac;

		}
		virtual ~HexadMesh()
		{
			release_mesh();
		}
	public:

		//Get element from container
		Hexahedron & handle_to_entity(HexadHandle    & handle);
		HalfFaceHH & handle_to_entity(HalfFaceHandle & handle);
		HalfEdgeHH & handle_to_entity(HalfEdgeHandle & handle);
		Vertex     & handle_to_entity(VertexHandle   & handle);
		Point      & handle_to_entity(PointHandle    & handle);

		//Next Prev Opposite Radial operator
		Hexahedron & next_hexahedron(HexadHandle   & handle);
		HalfFaceHH & next_half_face(HalfFaceHandle & handle);
		HalfEdgeHH & next_half_edge(HalfEdgeHandle & handle);
		Hexahedron & next_hexahedron(Hexahedron    & hh);
		HalfFaceHH & next_half_face(HalfFaceHH     & hf);
		HalfEdgeHH & next_half_edge(HalfEdgeHH     & he);


		Hexahedron & prev_hexahedron(HexadHandle   & handle);
		HalfFaceHH & prev_half_face(HalfFaceHandle & handle);
		HalfEdgeHH & prev_half_edge(HalfEdgeHandle & handle);
		Hexahedron & prev_hexahedron(Hexahedron    & hh);
		HalfFaceHH & prev_half_face(HalfFaceHH     & hf);
		HalfEdgeHH & prev_half_edge(HalfEdgeHH     & he);

		/// const

		//Get element from container
		const Hexahedron & handle_to_entity(HexadHandle    & handle) const;
		const HalfFaceHH & handle_to_entity(HalfFaceHandle & handle) const;
		const HalfEdgeHH & handle_to_entity(HalfEdgeHandle & handle) const;
		const Vertex     & handle_to_entity(VertexHandle   & handle) const;
		const Point      & handle_to_entity(PointHandle    & handle) const;

		//Next Prev Opposite Radial operator
		const Hexahedron & next_hexahedron(HexadHandle   & handle) const;
		const HalfFaceHH & next_half_face(HalfFaceHandle & handle) const;
		const HalfEdgeHH & next_half_edge(HalfEdgeHandle & handle) const;
		const Hexahedron & next_hexahedron(Hexahedron    & hh) const;
		const HalfFaceHH & next_half_face(HalfFaceHH     & hf) const;
		const HalfEdgeHH & next_half_edge(HalfEdgeHH     & he) const;


		const Hexahedron & prev_hexahedron(HexadHandle   & handle) const;
		const HalfFaceHH & prev_half_face(HalfFaceHandle & handle) const;
		const HalfEdgeHH & prev_half_edge(HalfEdgeHandle & handle) const;
		const Hexahedron & prev_hexahedron(Hexahedron    & hh) const;
		const HalfFaceHH & prev_half_face(HalfFaceHH     & hf) const;
		const HalfEdgeHH & prev_half_edge(HalfEdgeHH     & he) const;

		//Next Prev Opposite Radial operator
		HexadHandle     next_hexahedron_handle(HexadHandle   & handle);
		HalfFaceHandle  next_half_face_handle(HalfFaceHandle & handle);
		HalfEdgeHandle  next_half_edge_handle(HalfEdgeHandle & handle);
		HexadHandle     next_hexahedron_handle(Hexahedron    & hh);
		HalfFaceHandle  next_half_face_handle(HalfFaceHH     & hf);
		HalfEdgeHandle  next_half_edge_handle(HalfEdgeHH     & he);


		HexadHandle     prev_hexahedron_handle(HexadHandle   & handle);
		HalfFaceHandle  prev_half_face_handle(HalfFaceHandle & handle);
		HalfEdgeHandle  prev_half_edge_handle(HalfEdgeHandle & handle);
		HexadHandle     prev_hexahedron_handle(Hexahedron    & hh);
		HalfFaceHandle  prev_half_face_handle(HalfFaceHH     & hf);
		HalfEdgeHandle  prev_half_edge_handle(HalfEdgeHH     & he);

		HalfEdgeHH & radial_half_edge(HalfEdgeHandle   & handle);
		bool has_radial_half_edge(HalfEdgeHandle       & handle);
		HalfEdgeHH & mate_half_edge(HalfEdgeHandle     & handle);
		HalfFaceHH & opposite_half_face(HalfFaceHandle & handle);

		HalfEdgeHandle  radial_half_edge_handle(HalfEdgeHandle   & handle);
		HalfEdgeHandle  mate_half_edge_handle(HalfEdgeHandle     & handle);
		HalfFaceHandle  opposite_half_face_handle(HalfFaceHandle & handle);

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

		bool has_opposite_half_face(HalfFaceHandle& handle);
		bool has_vertex(HalfFaceHandle hfHandle, VertexHandle vHandle);
		bool has_point(HalfFaceHandle hfHandle, PointHandle pHandle);
		bool has_boundary_face(PolyhedronHandle & hh_);
		bool has_boundary_face(HedronIter & hIter_);

		/// is boundary
		bool is_boundary(HalfFaceHandle   & hfh_);
		bool is_boundary(VertexHandle     & vh_);
		bool is_boundary(PolyhedronHandle & hh_);
		bool is_boundary(HedronHandle & hh_, unsigned int type_);
		void update_boundary_point();

		Hexahedron & hexahedron(VertexHandle& handle);

		HalfEdgeHandle & half_edge_handle(HalfFaceHandle & fh, VertexHandle & vh);
		HalfFaceHandle & half_face_handle(HalfEdgeHandle & handle);

		/// Short-cut from halfedge to point
		PointHandle start_point_handle(HalfEdgeHandle& handle);
		Point& start_point(HalfEdgeHandle& handle);

		/// clear 
		void release_mesh();

		/// add hexahedron
		HexadHandle add_hedron(std::vector<Point>& pntVec);

		/// delete hexahedron
		HexadHandle erase_hedron(HexadHandle &handle);

		/// add point
		PointHandle add_point(Point &p);

		/// delete point
		int erase_point(PointHandle &handle);

		/// advanced operation
		int point_star(PointHandle ph, std::vector<HexadHandle> &hexadVec);
		int vertex_star(VertexHandle vh, std::vector<HexadHandle> &hexadVec);
		int edge_star(HalfEdgeHandle eh, std::vector<HexadHandle> &hexadVec);
		int vertex_half_edge_star(VertexHandle vh, std::vector<HalfEdgeHandle> &halfEdgeVec);
		int vertex_half_face_star(VertexHandle vh, std::vector<HalfFaceHandle> &halfFaceVec);

		/// build the topology
		bool build_topology();
		bool build_opposite_half_face();

		//-------------------------------------------------------------------------------------------------------------
		public:
			VertexHandle first_vertex_handle(HalfFaceHandle   fh_);
			VertexHandle first_vertex_handle(PolyhedronHandle hh_);
			VertexHandle first_vertex_handle(HedronIter    tIter_);
			VertexHandle next_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_);
			VertexHandle next_vertex_handle(HedronIter tIter_,    VertexHandle cvh_);
			VertexHandle prev_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_);
			VertexHandle prev_vertex_handle(HedronIter tIter_,    VertexHandle cvh_);

			HalfFaceHandle first_half_face_handle(HexadHandle    th_);
			HalfFaceHandle first_half_face_handle(HedronIter  tIter_);
			HalfEdgeHandle first_half_edge_handle(HalfFaceHandle fh_);

			VertexHandle from_vertex_handle(HalfEdgeHandle hh_);
			VertexHandle to_vertex_handle(HalfEdgeHandle   hh_);

			Hexahedron & iter_to_entity(HedronIter & iter);
			Vertex     & iter_to_entity(VertexIter & iter);

			void update_face_normals();
			Normal update_face_normal(HalfFaceHandle & hfh_);
			Normal normal(HalfFaceHandle & hfh_);


			/// access the point
			Point & iter_to_entity(PointIter  & iter);
			Point & point(PointHandle  & ph);
			Point & point(VertexHandle & handle);
			Point & point(PointIter    & iter);
			Point & point(VertexIter   & iter);

			/// element property
			int property(HalfFaceHH hf);
			int property(HalfFaceHandle hfh);
			
			void set_property(HalfFaceHH hf, int prop);
			void set_property(HalfFaceHandle hfh, int prop);
	//-----------------------------------------------------------------------------------------------------------------
	public:
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

				/** set the geometry position of a point handle
		*   \param in _ph the vertex handle
		*   \param in _p the position of the point
		*/
		//void set_point(const PointHandle & _ph, const Point & _p)
		//{
		//	_pc[_ph.idx()] = _p;
		//}

		void set_point(PointHandle _ph, Point _p)
		{
			_pc[_ph.idx()] = _p;
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
			return da_cos;
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

	//-----------------------------------------------------------------------------------------------------------------
	public:
		/** return the first iterator of the HalfFaceIter of the tetrahedron
		*  \param th_ the handle of the tetrahedron
		*  \return HedronFaceIter the iterator.
		*/
		HedronFaceIter hedron_face_iter(HexadHandle th_)
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

		/** return the first half edge handle of a hedron
		*   \param in hh_ the hedron handle
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
			return HedronIter(first_valid_hedron_handle().idx());
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

		/** return the hedron handle from the vertex
		*   \param vh the vertex handle
		*   \return PolyhedronHandle the handle of the Polyhedron
		*/
		PolyhedronHandle hedron_handle(VertexHandle vh)
		{
			if ((vh.idx() < 0) || (vh.idx() >= (int)_vc.size()))
			{
				return PolyhedronHandle(-1);
			}
			return HedronHandle(vh >> 3);
		}

		/** return the first adjacent hedron handle from the current hedron handle
		*   \param ph the polyhedron handle
		*   \return PolyhedronHandle the handle of the Polyhedron
		*/
		PolyhedronHandle first_adjacent_hedron_handle(PolyhedronHandle ph)
		{
			if ((ph.idx() < 0) || (ph.idx() >= (int)_hexac.size()))
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

		HalfEdgeHandle vertex_first_half_edge_handle(VertexHandle vh)
		{
			return HalfEdgeHandle((vh>>3)*24+vhecHH_[vh%8]);
		}
		//-----------------------------------------------------------------------------------------------------------------
		//API about _free_list

		bool is_valid(HexadHandle &hh_)
		{
			// valid handle and not in the free list
			if ( hh_>=0 && hh_<(int)_hexac.size() && find(_free_hedron_list.begin(), _free_hedron_list.end(), hh_) ==_free_hedron_list.end() )
				return true;
			else
				return false;
		}

		bool is_valid(HexadHandle &hh_) const
		{
			// valid handle and not in the free list
			if ( hh_>=0 && hh_<(int)_hexac.size() && find(_free_hedron_list.begin(), _free_hedron_list.end(), hh_) ==_free_hedron_list.end() )
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
			return is_valid(hedron_handle(vh_));
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
		HexadHandle next_valid_hedron_handle(HexadHandle &hh_)
		{
			//assert if input HedronHandle is out of range;
			assert ( (hh_.idx() < _npolyhedron ) && (hh_.idx() >= 0) );
			int idx = hh_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==_npolyhedron-1 )
					idx = 0;
				else 
					idx++;
				//find out valid hedron
				if( is_valid(HexadHandle(idx)) )
					return HexadHandle(idx);
			}
			while(idx!=hh_.idx());

			return HexadHandle(-1);
		}

		HexadHandle next_valid_hedron_handle(HexadHandle &hh_) const
		{
			//assert if input HedronHandle is out of range;
			assert ( (hh_.idx() < _npolyhedron ) && (hh_.idx() >= 0) );
			int idx = hh_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==_npolyhedron-1 )
					idx = 0;
				else 
					idx++;
				//find out valid HedronHandle
				if( is_valid(HexadHandle(idx)) )
					return HexadHandle(idx);
			}
			while(idx!=hh_.idx());

			return HexadHandle(-1);
		}


		//return prev valid handle
		//return -1 if no handle at all
		HexadHandle prev_valid_hedron_handle(HexadHandle &hh_)
		{
			//assert if input HedronHandle is out of range;
			assert ( (hh_.idx() < _npolyhedron ) && (hh_.idx() >= 0) );
			int idx = hh_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==0 )
					idx = _npolyhedron-1;
				else 
					idx--;
				//find out valid hedron
				if( is_valid(HexadHandle(idx)) )
					return HexadHandle(idx);
			}
			while(idx!=hh_.idx());

			return HexadHandle(-1);
		}

		HexadHandle prev_valid_hedron_handle(HexadHandle &hh_) const
		{
			//assert if input HedronHandle is out of range;
			assert ( (hh_.idx() < _npolyhedron ) && (hh_.idx() >= 0) );
			int idx = hh_.idx();
			//get rid of endless loop, can only loop once over the whole tetra handle
			do
			{
				if ( idx==0 )
					idx = _npolyhedron-1;
				else 
					idx--;
				//find out valid hedron
				if( is_valid(HexadHandle(idx)) )
					return HexadHandle(idx);
			}
			while(idx!=hh_.idx());

			return HexadHandle(-1);
		}

		//return first valid tetrahedron handle
		HexadHandle first_valid_hedron_handle()
		{
			return next_valid_hedron_handle(HexadHandle(_npolyhedron-1));
		}

		HexadHandle last_valid_tetrahedron_handle()
		{
			return prev_valid_hedron_handle(HexadHandle(0));
		}

		void add_to_free_hedron_list(HexadHandle &hh_)
		{
			if (find(_free_hedron_list.begin(), _free_hedron_list.end(), hh_) == _free_hedron_list.end())
				_free_hedron_list.push_back(hh_);
		}

		void add_to_free_point_list(PointHandle &ph_)
		{
			if (find(_free_point_list.begin(), _free_point_list.end(), ph_) == _free_point_list.end())
				_free_point_list.push_back(ph_);
			//#需要维护pvm
		}


		//Get first free hedron location, or return -1 if empty
		HexadHandle pop_from_free_tetra_list()
		{
			if( _free_hedron_list.size()>0 )
				return *_free_hedron_list.begin();
			else
				return HexadHandle(-1);
		}

		//Get first free hedron location, then erase it. or return -1 if empty
		HexadHandle pop_erase_from_free_tetra_list()
		{
			if ( _free_hedron_list.size()>0 )
			{
				HexadHandle hh = *_free_hedron_list.begin();
				_free_hedron_list.erase(_free_hedron_list.begin());
				return hh;
			}
			else
				return HexadHandle(-1);
		}


		PointHandle pop_from_free_point_list()
		{
			if( _free_point_list.size()>0 )
				return *_free_point_list.begin();
			else
				return PointHandle(-1);
			//#需要维护pvm
		}

		PointHandle pop_erase_from_free_point_list()
		{
			if ( _free_point_list.size()>0 )
			{
				PointHandle ph = *_free_point_list.begin();
				_free_point_list.erase(_free_point_list.begin());
				_pvm[ph].clear();
				return ph;
			}
			else
				return PointHandle(-1);
			//#需要维护pvm
		}

		// find whether there is a valid point in point vector 
		PointHandle find_point(PointFindIF<Point>& pif)
		{
			std::vector<Point>::iterator iter=_pc.begin();
			while( (iter = find_if(iter,_pc.end(),pif)) != _pc.end() )
			{
				if ( is_valid(PointHandle(iter-_pc.begin())) )
					return PointHandle(iter-_pc.begin());
				else
				{
					iter++;
					if ( iter==_pc.end() )
						break;
				}
			}
			return PointHandle(-1);
		}

		// valid hedron number
		int size_valid_hedron()
		{
			return _npolyhedron - (int)_free_hedron_list.size();
		}
		// clean _free_hedron_list and _free_point_list & filling the vacant
		void clean_garbage();

	//-----------------------------------------------------------------------------------------------------------------
	/** 
	* members of the HexadMesh
	*/
	protected:
		HexadContainer _hexac;
		HFContainer _hfc;
		HEContainer _hec;
		std::vector<HedronHandle> _free_hedron_list;
		std::vector<PointHandle> _free_point_list;
	};

//---------------------------------------------------------------------------------------------------------------------
} // namespace VolumMesh
//---------------------------------------------------------------------------------------------------------------------
#endif 