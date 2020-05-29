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
* this file defines the circulators                                                                                   *
* $ written by xianc $                                                                                                *
* $ Version 0.1 $                                                                                                     *
* $ Date: 2010-02-08 $                                                                                                *
**********************************************************************************************************************/
/*--------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		      *
*					   Email: chuhuaxian@gmail.com																      *
*--------------------------------------------------------------------------------------------------------------------*/ 
/*---------------------------------------------------------------------------------------------------------------------
* Modified Information                                                                                                *
* Modified by Chuhua Xian                                                                                             *
* Email : chuhuaxian@gmail.com                                                                                        *
* Modified Date : 2010.8.7.                                                                                           *
*--------------------------------------------------------------------------------------------------------------------*/
#ifndef _VOLUME_MESH_CIRCULATOR_
#define _VOLUME_MESH_CIRCULATOR_

#include <VolumeMesh/system/config.h>
#include <VolumeMesh/Mesh/Handles.h>

#include <set>
#include <assert.h>

namespace VolumeMesh
{
	namespace Circulator
	{
		/** A HedronHedronIter is used to circulate over all polyhedrons of a polyhedron,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class HedronHedronIterT
		{
		public:
			/// default constructor
			HedronHedronIterT() : mesh_(NULL), lap_counter_(0)
			{
			}

			/// construct with mesh and half face handle
			HedronHedronIterT(MeshT & _mesh, HalfFaceHandle _start, bool _end = false) : mesh_(&_mesh), 
				              start_(_mesh.hedron_handle(mesh_->opposite_half_face_handle(_start))), 
							  ph_(start_), hfh_(_start), lap_counter_(_end)
			{
				if (!ph_.is_valid())
				{
					do
					{
						hfh_ = mesh_->next_half_face_handle(hfh_);
						ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
					}while((!ph_.is_valid()) && ((_start != hfh_)));
					start_ = ph_;
				}
			}

			/// construct with mesh and polyhedron handle
			HedronHedronIterT(MeshT & _mesh, PolyhedronHandle _ph, bool _end = false) : mesh_ (&_mesh),
				              start_(_mesh.first_adjacent_hedron_handle(_ph)),
				              ph_(start_),
							  hfh_(_mesh.first_half_face_handle(_ph)),
							  lap_counter_(_end)
			{
				if (!ph_.is_valid())
				{
					do
					{
						hfh_ = mesh_->next_half_face_handle(hfh_);
						ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
					}while((!ph_.is_valid()) && ((_start != hfh_)));
					start_ = ph_;
				}
			}

			/// copy constructor
			HedronHedronIterT(const HedronHedronIterT & _rhs) : mesh_(_rhs.mesh_), 
				start_(_rhs.start_), 
				ph_(_rhs.ph_),
				hfh_(_rhs.hfh_),
				lap_counter_(_rhs.lap_counter_)
			{
			}

			~ HedronHedronIterT()
			{
			}
		public:
			HedronHedronIterT & operator = (const HedronHedronIterT<MeshT> & _rhs)
			{
				mesh_  = _rhs.mesh_;
				start_ = _rhs.start_;
				ph_    = _rhs.ph_;
				hfh_   = _rhs.hfh_;

				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const HedronHedronIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)       &&
					   (start_       == _rhs.start_)      &&
					   (ph_          == _rhs.ph_)         &&
					   (hfh_         == _rhs.hfh_)        &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const HedronHedronIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			HedronHedronIterT & operator ++ () 
			{
				assert(mesh_); 
				hfh_ = mesh_->next_half_face_handle(hfh_);
				ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
				while((!ph_.is_valid()) && ((start_ != ph_) || (lap_counter_ == 0)))
				{
					hfh_ = mesh_->next_half_face_handle(hfh_);
					ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
				}

				if (start_ == ph_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			HedronHedronIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == ph_)
				{
					-- lap_counter_;
				} 
				hfh_ = mesh_->prev_half_face_handle(hfh_);
				ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
				while((!ph_.is_valid()) && ((start_ != ph_) || (lap_counter_ == 0))) 
				{
					//if (start_ == ph_)
					//{
					//	-- lap_counter_;
					//} 
					hfh_ = mesh_->prev_half_face_handle(hfh_);
					ph_  = mesh_->hedron_handle(mesh_->opposite_half_face_handle(hfh_));
				}
				return *this;
			}

			/// return current half face handle
			PolyhedronHandle & current_hedron_handle()
			{
				return ph_;
			}

			/// return the polyhedron handle
			PolyhedronHandle handle() const
			{
				//assert(mesh_);
				return ph_;
			}

			/// return a reference HalfFace 
			Polyhedron & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(ph_);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return ph_.is_valid() && ((start_ != ph_) || (lap_counter_ == 0));
			}

		protected:
			MeshT             * mesh_;
			PolyhedronHandle    start_, ph_;
			HalfFaceHandle      hfh_;
			int                 lap_counter_;
		};

		//-------------------------------------------------------------------------------------------------------------

		/** A HedronFaceIter is used to circulate over all half face of a polyhedron,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class HedronHalfFaceIterT
		{
		public:
			/// default constructor
			HedronHalfFaceIterT() : mesh_(NULL), lap_counter_(0)
			{

			}

			/// construct with mesh and polyhedron
			HedronHalfFaceIterT(MeshT & _mesh, PolyhedronHandle _start, bool _end = false) : mesh_(&_mesh), 
				            start_(_mesh.first_half_face_handle(_start)), hfh_(start_), lap_counter_(_end)
			{
			}

			/// construct with mesh and half face
			HedronHalfFaceIterT(MeshT & _mesh, HalfFaceHandle _hfh, bool _end = false) : mesh_ (&_mesh),
				            start_(_hfh),
							hfh_(_hfh),
							lap_counter_(_end)
			{
			}

			/// copy constructor
			HedronHalfFaceIterT(const HedronHalfFaceIterT & _rhs) : mesh_(_rhs.mesh_), 
				                                           start_(_rhs.start_), 
														   hfh_(_rhs.hfh_),
														   lap_counter_(_rhs.lap_counter_)
			{
			}

			~ HedronHalfFaceIterT()
			{
			}
		public:
			HedronHalfFaceIterT & operator = (const HedronHalfFaceIterT<MeshT> & _rhs)
			{
				mesh_  = _rhs.mesh_;
				start_ = _rhs.start_;
				hfh_   = _rhs.hfh_;
				
				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const HedronHalfFaceIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)       &&
					   (start_       == _rhs.start_)      &&
					   (hfh_         == _rhs.hfh_)        &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const HedronHalfFaceIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			HedronHalfFaceIterT & operator ++ () 
			{
				assert(mesh_); 
				hfh_ = mesh_->next_half_face_handle(hfh_);
				if (start_ == hfh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			HedronHalfFaceIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == hfh_)
				{
					-- lap_counter_;
				} 
				hfh_ = mesh_->prev_half_face_handle(hfh_);
				return *this;
			}

			/// return current half face handle
			HalfFaceHandle & current_half_face_handle()
			{
				return hfh_;
			}

			/// return the polyhedron handle
			HalfFaceHandle handle() const
			{
				//assert(mesh_);
				return hfh_;
			}

			/// return a reference HalfFace 
			HalfFace & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(hfh_);
			}

			/** 
			 * return if the circulator is still valid, and after the circulator round around , it becoms invalid
			 */
			operator bool() const
			{
				return hfh_.is_valid() && ((start_ != hfh_) || (lap_counter_ == 0));
			}

		protected:
			MeshT          * mesh_;
			HalfFaceHandle   start_, hfh_;			
			int              lap_counter_;
		};

		//-------------------------------------------------------------------------------------------------------------

		/** A HedronFaceIter is used to circulate over all half face of a polyhedron,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class HedronHalfEdgeIterT
		{
		public:
			/// default constructor
			HedronHalfEdgeIterT() : mesh_(NULL), lap_counter_(0)
			{

			}

			/// construct with mesh and polyhedron
			HedronHalfEdgeIterT(MeshT & _mesh, typename MeshT::HedronHandle _ph, bool _end = false) : mesh_(&_mesh), ph_(_ph),
				                start_(_mesh.first_hedron_half_edge_handle(_ph)), heh_(start_), lap_counter_(_end)
			{
			}

			/// construct with mesh and half face


			/// copy constructor
			HedronHalfEdgeIterT(const HedronHalfEdgeIterT & _rhs) : mesh_(_rhs.mesh_), ph_(_rhs.ph_),
				start_(_rhs.start_), 
				heh_(_rhs.heh_),
				lap_counter_(_rhs.lap_counter_)
			{
			}

			~ HedronHalfEdgeIterT()
			{
			}
		public:
			HedronHalfEdgeIterT & operator = (const HedronHalfEdgeIterT & _rhs)
			{
				mesh_  = _rhs.mesh_;
				ph_    = _rhs.ph_;
				start_ = _rhs.start_;				
				heh_   = _rhs.heh_;

				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const HedronHalfEdgeIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)       &&
					   (ph_          == _rhs.ph_)         &&
					   (start_       == _rhs.start_)      &&
					   (heh_         == _rhs.heh_)        &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const HedronHalfEdgeIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			HedronHalfEdgeIterT & operator ++ () 
			{
				assert(mesh_); 
				heh_ = mesh_->next_hedron_half_edge_handle(ph_, heh_);
				if (start_ == heh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			HedronHalfEdgeIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == heh_)
				{
					-- lap_counter_;
				} 
				heh_ = mesh_->prev_hedron_half_edge_handle(ph_, heh_);
				return *this;
			}

			/// return current half face handle
			HalfEdgeHandle & current_half_edge_handle()
			{
				return heh_;
			}

			/// return the polyhedron handle
			HalfEdgeHandle handle() const
			{
				//assert(mesh_);
				return heh_;
			}

			/// return a reference HalfEdge 
			HalfEdge & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(heh_);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becomes invalid
			*/
			operator bool() const
			{
				return heh_.is_valid() && ((start_ != heh_) || (lap_counter_ == 0));
			}

		protected:
			MeshT                         * mesh_;
			typename MeshT::HedronHandle    ph_;
			HalfEdgeHandle                  start_, heh_;			
			int                             lap_counter_;
		};


		//-----------------------------------------------------------------------------------------------------------------

		/** A HalfFaceVertexIter is used to circulate over all vertex of a half face,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class HalfFaceVertexIterT
		{
		public:
			/// default constructor
			HalfFaceVertexIterT() : mesh_(NULL), lap_counter_(0)
			{

			}

			/// construct with mesh and polyhedron
			HalfFaceVertexIterT(MeshT & _mesh, HalfFaceHandle _start, bool _end = false) : mesh_(&_mesh), 
				start_(_mesh.first_half_edge_handle(_start)), lap_counter_(_end)
			{
				heh_ = start_;
			}

			/// construct with mesh and half face
			HalfFaceVertexIterT(MeshT & _mesh, HalfEdgeHandle _heh, bool _end = false) : mesh_ (&_mesh),
				                start_(_heh),
				                heh_(_heh),
								lap_counter_(_end)
			{
			}

			/// copy constructor
			HalfFaceVertexIterT(const HalfFaceVertexIterT & _rhs) : mesh_(_rhs.mesh_), 
				                start_(_rhs.start_),
								heh_(_rhs.heh_),
								lap_counter_(_rhs.lap_counter_)
			{
			}

			~ HalfFaceVertexIterT()
			{
			}
		public:
			HalfFaceVertexIterT & operator = (const HalfFaceVertexIterT<MeshT> _rhs)
			{
				mesh_  = _rhs.mesh_;
				start_ = _rhs.start_;
				heh_   = _rhs.heh_;

				lap_counter_ = _rhs.lap_counter_;

				return *this;
			}


			bool operator == (const HalfFaceVertexIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)  &&
					   (start_       == _rhs.start_) &&
					   (heh_         == _rhs.heh_)   &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const HalfFaceVertexIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			HalfFaceVertexIterT & operator ++ () 
			{
				assert(mesh_); 
				heh_ = mesh_->next_half_edge_handle(heh_);
				if (start_ == heh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			HalfFaceVertexIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == heh_)
				{
					-- lap_counter_;
				} 
				heh_ = mesh_->prev_half_face_handle(heh_);
				return *this;
			}

			/// return current half face handle
			HalfEdgeHandle & current_half_edge_handle()
			{
				return heh_;
			}

			/// return the polyhedron handle
			VertexHandle handle() const
			{
				assert(mesh_);
				return mesh_->from_vertex_handle(heh_);
			}

			/// return a reference vertex 
			Vertex & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(handle());
			}

			/// return a vertex pointer
			Vertex * operator-> () const
			{
				assert(mesh_);
				return &mesh_->handle_to_entity(handle());
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return heh_.is_valid() && ((start_ != heh_) || (lap_counter_ == 0));
			}

		protected:
			MeshT          * mesh_;
			HalfEdgeHandle   start_, heh_;			
			int              lap_counter_;
		};
		//-----------------------------------------------------------------------------------------------------------------

		/** A FaceHalfedgeIter is used to circulate over all halfedges of a half face,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class FaceHalfedgeIterT
		{
		public:
			/// default constructor
			FaceHalfedgeIterT() : mesh_(NULL), lap_counter_(0)
			{
			}

			/// construct with mesh and polyhedron
			FaceHalfedgeIterT(MeshT & _mesh, HalfFaceHandle _start, bool _end = false) : mesh_(&_mesh), 
				start_(_mesh.first_half_edge_handle(_start)), lap_counter_(_end)
			{
				heh_ = start_;
			}

			/// construct with mesh and half face
			FaceHalfedgeIterT(MeshT & _mesh, HalfEdgeHandle _heh, bool _end = false) : mesh_ (&_mesh),
				start_(_heh),
				heh_(_heh),
				lap_counter_(_end)
			{
			}

			/// copy constructor
			FaceHalfedgeIterT(const FaceHalfedgeIterT & _rhs) : mesh_(_rhs.mesh_), 
				start_(_rhs.start_),
				heh_(_rhs.heh_),
				lap_counter_(_rhs.lap_counter_)
			{
			}

			~ FaceHalfedgeIterT()
			{
			}

		public:
			FaceHalfedgeIterT & operator = (const FaceHalfedgeIterT<MeshT> _rhs)
			{
				mesh_  = _rhs.mesh_;
				start_ = _rhs.start_;
				heh_   = _rhs.heh_;

				lap_counter_ = _rhs.lap_counter_;

				return *this;
			}


			bool operator == (const FaceHalfedgeIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)  &&
					   (start_       == _rhs.start_) &&
					   (heh_         == _rhs.heh_)   &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const FaceHalfedgeIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			FaceHalfedgeIterT & operator ++ () 
			{
				assert(mesh_); 
				heh_ = mesh_->next_half_edge_handle(heh_);
				if (start_ == heh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			FaceHalfedgeIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == heh_)
				{
					-- lap_counter_;
				} 
				heh_ = mesh_->prev_half_face_handle(heh_);
				return *this;
			}

			/// return current half face handle
			FaceHalfedgeIterT & current_half_edge_handle()
			{
				return heh_;
			}

			/// return the polyhedron handle
			HalfEdgeHandle handle() const
			{
				assert(mesh_);
				return heh_;
			}

			/// return a reference vertex 
			HalfEdge & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(handle());
			}

			/// return a vertex pointer
			HalfEdge * operator-> () const
			{
				assert(mesh_);
				return &mesh_->handle_to_entity(handle());
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return heh_.is_valid() && ((start_ != heh_) || (lap_counter_ == 0));
			}

		protected:
			MeshT          * mesh_;
			HalfEdgeHandle   start_, heh_;			
			int              lap_counter_;
		};


		//-------------------------------------------------------------------------------------------------------------

		/** A HedronVertexIter is used to circulate over all vertex of a polyhedron,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class HedronVertexIterT
		{
		public:
			/// default constructor
			HedronVertexIterT() : mesh_(NULL), lap_counter_(0)
			{

			}

			/// construct with mesh and polyhedron
			HedronVertexIterT(MeshT & _mesh, PolyhedronHandle _hh, HalfFaceHandle _start, bool _end = false) : mesh_(&_mesh), 
				             hh_(_hh),							 
							 start_(_mesh.first_vertex_handle(mesh_->first_half_edge_handle(_start))), 
							 lap_counter_(_end)
			{
				vh_ = start_;
			}

			/// construct with mesh and half face
			HedronVertexIterT(MeshT & _mesh, PolyhedronHandle _hh, VertexHandle _vh, bool _end = false) : mesh_ (&_mesh),
				start_(_vh),
				vh_(_vh),
				hh_(_hh),
				lap_counter_(_end)
			{
			}

			/// copy constructor
			HedronVertexIterT(const HedronVertexIterT & _rhs) : mesh_(_rhs.mesh_), 
				              start_(_rhs.start_),
				              vh_(_rhs.vh_),
				              hh_(_rhs.hh_),
				              lap_counter_(_rhs.lap_counter_)
			{
			}

			~ HedronVertexIterT()
			{
			}
		public:
			HedronVertexIterT & operator = (const HedronVertexIterT<MeshT> _rhs)
			{
				mesh_  = _rhs.mesh_;
				start_ = _rhs.start_;
				vh_    = _rhs.vh_;
				hh_    = _rhs.hh_;

				lap_counter_ = _rhs.lap_counter_;

				return *this;
			}


			bool operator == (const HedronVertexIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)  &&
					   (start_       == _rhs.start_) &&
					   (vh_          == _rhs.vh_)   &&
					   (hh_          == _rhs.hh_)   &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const HedronVertexIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			HedronVertexIterT & operator ++ () 
			{
				assert(mesh_); 
				vh_ = mesh_->next_vertex_handle(hh_, vh_);
				if (start_ == vh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			HedronVertexIterT & operator -- () 
			{
				assert(mesh_);
				if (start_ == vh_)
				{
					-- lap_counter_;
				} 
				vh_ = mesh_->prev_half_face_handle(hh_, vh_);
				return *this;
			}

			/// return current vertex handle
			VertexHandle & current_vertex_handle()
			{
				return vh_;
			}

			/// return the vertex handle
			VertexHandle handle() const
			{
				assert(mesh_);
				return vh_;
			}

			/// return a reference vertex 
			Vertex & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(handle());
			}

			/// return a vertex pointer
			Vertex * operator-> () const
			{
				assert(mesh_);
				return &mesh_->handle_to_entity(handle());
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return vh_.is_valid() && ((start_ != vh_) || (lap_counter_ == 0));
			}

		protected:
			MeshT          * mesh_;
			VertexHandle     start_, vh_;
			PolyhedronHandle hh_;
			int              lap_counter_;
		};

		//-------------------------------------------------------------------------------------------------------------

		/** A VertexHedronIter is used to circulate over all polyhedrons of a Vertex,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class VertexHedronIterT
		{
		public:
			/// default constructor
			VertexHedronIterT() : mesh_(NULL), lap_counter_(0)
			{
			}

			/// construct with mesh and polyhedron handle
			VertexHedronIterT(MeshT & _mesh, const std::vector<PolyhedronHandle> & _hedrons, bool _end = false) :
			                  mesh_ (&_mesh),
				              hedrons_(_hedrons),
							  ph_iter_(0),
							  lap_counter_(_end)
			{
			}

			/// copy constructor
			VertexHedronIterT(const VertexHedronIterT & _rhs) : mesh_(_rhs.mesh_), 
				                                                ph_iter_(_rhs.ph_iter_), 
																hedrons_(_rhs.hedrons_),
																lap_counter_(_rhs.lap_counter_)
			{
			}

			~ VertexHedronIterT()
			{
			}
		public:
			VertexHedronIterT & operator = (const VertexHedronIterT<MeshT> _rhs)
			{
				mesh_        = _rhs.mesh_;
				hedrons_     = _rhs.hedrons_;
				ph_iter_     = _rhs.ph_iter_;
				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const VertexHedronIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)        &&
					   (hedrons_     == _rhs.hedrons_)     &&
					   (ph_iter_     == _rhs.ph_iter_)     &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const VertexHedronIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			VertexHedronIterT & operator ++ () 
			{
				assert(mesh_);
				++ ph_iter_;
				if (ph_iter_ == 0)
				{
					++ lap_counter_;
				}

				return *this;
			}

			/// decrement operator
			VertexHedronIterT & operator -- () 
			{
				assert(mesh_);
				if (ph_iter_ == hedrons_.size())
				{
					-- lap_counter_;
				}
				-- ph_iter_;
				return *this;
			}

			/// return current half face handle
			PolyhedronHandle & current_hedron_handle()
			{
				return hedrons_[ph_iter_];
			}

			/// return the polyhedron handle
			PolyhedronHandle handle() const
			{
				//assert(mesh_);
				return hedrons_[ph_iter_];
			}

			/// return a reference HalfFace 
			Polyhedron & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(hedrons_[ph_iter_]);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return (ph_iter_ != hedrons_.size()) && ((ph_iter_ != 0) || (lap_counter_ == 0));
			}

		protected:
			MeshT                                  * mesh_;
			std::vector<PolyhedronHandle>            hedrons_;
			int                                      ph_iter_;
			int                                      lap_counter_;
		};

		//-------------------------------------------------------------------------------------------------------------
		
		/** A VertexHalfFaceIter is used to circulate over all HalfFaces of a Vertex,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class VertexHalfFaceIterT
		{
		public:
			VertexHalfFaceIterT() : mesh_(NULL), lap_counter_(0)
			{
			}
			VertexHalfFaceIterT(MeshT & _mesh, VertexHandle _vh, HalfFaceHandle _hfh, bool _end = false) : 
			    mesh_(&_mesh),
				vh_(_vh),
				start_(_hfh),
				hfh_(_hfh),
				lap_counter_(_end)
			{
			}
			VertexHalfFaceIterT(const VertexHalfFaceIterT & _rhs) : 
				mesh_(_rhs.mesh_), 
				vh_(_rhs.vh_), 
				start_(_rhs.start_),
				hfh_(_rhs.hfh_),
				lap_counter_(_rhs.lap_counter_)
			{
			}
			~VertexHalfFaceIterT(){}

		public:
			VertexHalfFaceIterT & operator = (const VertexHalfFaceIterT<MeshT> _rhs)
			{
				mesh_        = _rhs.mesh_;
				vh_          = _rhs.vh_;
				start_       = _rhs.start_;
				hfh_         = _rhs.hfh_;
				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const VertexHalfFaceIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)       &&
					   (vh_          == _rhs.vh_)         &&
					   (start_         == _rhs.start_)    &&
					   (hfh_         == _rhs.hfh_)        &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const VertexHalfFaceIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			VertexHalfFaceIterT & operator ++ () 
			{
				assert(mesh_);
				hfh_ = mesh_->next_vertex_half_face_handle(vh_, hfh_);
				if (start_ == hfh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			VertexHalfFaceIterT & operator -- () 
			{
				assert(mesh_);
				hfh_ = mesh_->prev_vertex_half_face_handle(vh_, hfh_);
				if (start_ == hfh_)
				{
					-- lap_counter_;
				}
				return *this;
			}

			/// return current half face handle
			HalfFaceHandle & current_half_face_handle()
			{
				return hfh_;
			}

			/// return the polyhedron handle
			HalfFaceHandle handle() const
			{
				//assert(mesh_);
				return hfh_;
			}

			/// return a reference HalfFace 
			HalfFace & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(hfh_);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return (start_ != hfh_) || (lap_counter_ == 0);
			}

		protected:
			MeshT            * mesh_;
			VertexHandle       vh_;
			HalfFaceHandle     start_, hfh_;
			int                lap_counter_;
		};

		
		//-------------------------------------------------------------------------------------------------------------
		
		/** A VertexHalfEdgeIter is used to circulate over all HalfEdges of a Vertex,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class VertexHalfEdgeIterT
		{
		public:
			VertexHalfEdgeIterT() : mesh_(NULL), loc_flag_(0), lap_counter_(0)
			{
			}
			VertexHalfEdgeIterT(MeshT & _mesh, VertexHandle _vh, HalfEdgeHandle _heh, int _loc_flag = 0, bool _end = false) : 
			    mesh_(&_mesh),
				vh_(_vh),
				start_(_heh),
				heh_(_heh),
				loc_flag_(_loc_flag),
				lap_counter_(_end)
			{
			}
			VertexHalfEdgeIterT(const VertexHalfEdgeIterT & _rhs) : 
				mesh_(_rhs.mesh_), 
				vh_(_rhs.vh_), 
				start_(_rhs.start_),
				heh_(_rhs.heh_),
				loc_flag_(_rhs.loc_flag_),
				lap_counter_(_rhs.lap_counter_)
			{
			}
			~VertexHalfEdgeIterT(){}

		public:
			VertexHalfEdgeIterT & operator = (const VertexHalfEdgeIterT<MeshT> _rhs)
			{
				mesh_        = _rhs.mesh_;
				vh_          = _rhs.vh_;
				start_       = _rhs.start_;
				heh_         = _rhs.heh_;
				loc_flag_    = _rhs.loc_flag_;
				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const VertexHalfEdgeIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)        &&
					   (vh_          == _rhs.vh_)          &&
					   (start_       == _rhs.start_)       &&
					   (heh_         == _rhs.heh_)         &&
					   (loc_flag_    == _rhs.loc_flag_)    &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const VertexHalfEdgeIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			VertexHalfEdgeIterT & operator ++ () 
			{
				assert(mesh_);
				if (loc_flag_)
				{
					heh_ = mesh_->next_half_edge_handle(heh_);
					loc_flag_ = 0;
				}
				else
				{
					heh_ = mesh_->mate_half_edge_handle(heh_);
					loc_flag_ = 1;
				}
				//!loc_flag_;
				if (start_ == heh_)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			VertexHalfEdgeIterT & operator -- () 
			{
				assert(mesh_);
				if (loc_flag_)
				{
					heh_ = mesh_->mate_half_edge_handle(heh_);
					loc_flag_ = 0;
				}
				else
				{
					heh_ = mesh_->prev_half_edge_handle(heh_);
					loc_flag_ = 1;
				}
				//!loc_flag_;
				if (start_ == heh_)
				{
					-- lap_counter_;
				}
				return *this;
			}

			/// return current half edge handle
			HalfEdgeHandle & current_half_edge_handle()
			{
				return heh_;
			}

			/// return the half edge handle
			HalfEdgeHandle handle() const
			{
				//assert(mesh_);
				return heh_;
			}

			/// return a reference HalfEdge 
			HalfFace & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(heh_);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return (heh_ != start_) || (lap_counter_ == 0);
			}

		protected:
			MeshT            * mesh_;
			VertexHandle       vh_;
			HalfEdgeHandle     start_, heh_;
			int                loc_flag_;
			int                lap_counter_;
		};

		//-------------------------------------------------------------------------------------------------------------

		/** A PointHedronIter is used to circulate over all topological hedrons of a point,
		just like STL iterator.  Its
		operator* returns the current corresponding entity.
		*/
		template <class MeshT>
		class PointHedronIterT
		{
		public:
			/// default constructor
			PointHedronIterT() : mesh_(NULL), lap_counter_(0)
			{
			}

			/// construct with mesh and polyhedron handle
			PointHedronIterT(MeshT & _mesh, const std::vector<PolyhedronHandle> & _hedrons, bool _end = false) : 
			                 mesh_ (&_mesh),
				             hedrons_(_hedrons),
							 ph_iter_(0),
							 lap_counter_(_end)
			{
			}

			/// copy constructor
			PointHedronIterT(const PointHedronIterT & _rhs) : mesh_(_rhs.mesh_), 
				                                              hedrons_(_rhs.hedrons_), 
															  ph_iter_(_rhs.ph_iter_),
															  lap_counter_(_rhs.lap_counter_)
			{
			}

			~ PointHedronIterT()
			{
			}
		public:
			PointHedronIterT & operator = (const PointHedronIterT<MeshT> _rhs)
			{
				mesh_        = _rhs.mesh_;
				hedrons_     = _rhs.hedrons_;
				ph_iter_     = _rhs.ph_iter_;
				lap_counter_ = _rhs.lap_counter_;
				return *this;
			}
			bool operator == (const PointHedronIterT & _rhs) const
			{
				return (mesh_        == _rhs.mesh_)        &&
					   (hedrons_     == _rhs.hedrons_)     &&
					   (ph_iter_     == _rhs.ph_iter_)     &&
					   (lap_counter_ == _rhs.lap_counter_);
			}
			bool operator != (const PointHedronIterT & _rhs) const
			{
				return !operator == (_rhs);
			}

			/// increment operator
			PointHedronIterT & operator ++ () 
			{
				assert(mesh_);
				++ ph_iter_;
				if (ph_iter_ == 0)
				{
					++ lap_counter_;
				}
				return *this;
			}

			/// decrement operator
			PointHedronIterT & operator -- () 
			{
				assert(mesh_);
				if (ph_iter_ == hedrons_.size())
				{
					-- lap_counter_;
				}
				-- ph_iter_;
				return *this;
			}

			/// return current half face handle
			PolyhedronHandle & current_hedron_handle()
			{
				return hedrons_[ph_iter_];
			}

			/// return the polyhedron handle
			PolyhedronHandle handle() const
			{
				//assert(mesh_);
				return hedrons_[ph_iter_];
			}

			/// return a reference HalfFace 
			Polyhedron & operator * () const
			{
				assert(mesh_);
				return mesh_->handle_to_entity(hedrons_[ph_iter_]);
			}

			/** 
			* return if the circulator is still valid, and after the circulator round around , it becoms invalid
			*/
			operator bool() const
			{
				return (ph_iter_ != hedrons_.size()) && ((ph_iter_ != 0) || (lap_counter_ == 0));
			}

		protected:
			MeshT                                  * mesh_;
			std::vector<PolyhedronHandle>            hedrons_;
			int                                      ph_iter_;
			int                                      lap_counter_;
		};





	}

}

//---------------------------------------------------------------------------------------------------------------------

#endif