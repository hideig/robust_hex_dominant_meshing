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


#include <VolumeMesh/Mesh/TetraMesh.h>
#include <algorithm>
#include <set>

namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	TetraMesh::TetraMesh() : BaseMesh(0x0001)
	{

	}

	TetraMesh::~TetraMesh()
	{
		release_mesh();
	}

	/** Retrieve Tetra stored in the TetraMesh, translate the handle to the item
	*	\param  handle the index of the Tetra to retrieve, must be valid. 
	*   \return Tetrahedron which indicated by the handle
	*/
	Tetrahedron & TetraMesh::handle_to_entity(TetraHandle &handle)
	{
		if ( ((unsigned int)handle.idx() < _tetc.size()) && (handle.idx() >= 0))
		{
			return _tetc[handle.idx()];
		}
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve Tetra stored in the TetraMesh, translate the handle to the item
	*	\param  handle the index of the Tetra to retrieve, must be valid. 
	*   \return Tetrahedron which indicated by the handle
	*/
	const Tetrahedron & TetraMesh::handle_to_entity(TetraHandle & handle) const
	{
		if ( ((unsigned int)handle.idx() < _tetc.size()) && (handle.idx() >= 0))
		{
			return _tetc[handle.idx()];
		}
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH to retrieve, must be valid. 
	*   \return HalfFaceTH which indicated by the handle
	*/
	HalfFaceTH& TetraMesh::handle_to_entity(HalfFaceHandle &handle)
	{
		if (((unsigned int)handle.idx() < _hfc.size()) && (handle.idx() >= 0))
		{
			return _hfc[handle.idx()];
		} 
		else
		{
			return *((HalfFaceTH*)NULL);
		}
	}

	/** Retrieve HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH to retrieve, must be valid. 
	*   \return HalfFaceTH which indicated by the handle
	*/
	const HalfFaceTH& TetraMesh::handle_to_entity(HalfFaceHandle &handle) const
	{
		if (((unsigned int)handle.idx() < _hfc.size()) && (handle.idx() >= 0))
		{
			return _hfc[handle.idx()];
		} 
		else
		{
			return *((HalfFaceTH*)NULL);
		}
	}

	/** Retrieve HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH to retrieve, must be valid. 
	*   \return HalfEdgeTH which indicated by the handle
	*/
	HalfEdgeTH& TetraMesh::handle_to_entity(HalfEdgeHandle &handle)
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0))
		{
			return _hec[handle.idx()];
		} 
		else
		{
			return *((HalfEdgeTH*)NULL);
		}
	}

	/** Retrieve HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH to retrieve, must be valid. 
	*   \return HalfEdgeTH which indicated by the handle
	*/
	const HalfEdgeTH& TetraMesh::handle_to_entity(HalfEdgeHandle &handle) const
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0))
		{
			return _hec[handle.idx()];
		} 
		else
		{
			return *((HalfEdgeTH*)NULL);
		}
	}

	/** Retrieve Vertex stored in the TetraMesh.
	*	\param  handle the index of the Vertex to retrieve, must be valid. 
	*   \return Vertex which indicated by the handle
	*/
	Vertex& TetraMesh::handle_to_entity(VertexHandle& handle)
	{
		if (((unsigned int)handle.idx() < _vc.size()) && (handle.idx() >= 0))
		{
			return _vc[handle.idx()];
		} 
		else
		{
			return *((Vertex*)NULL);
		}
	}

	/** Retrieve Vertex stored in the TetraMesh.
	*	\param  handle the index of the Vertex to retrieve, must be valid. 
	*   \return Vertex which indicated by the handle
	*/
	const Vertex& TetraMesh::handle_to_entity(VertexHandle& handle) const
	{
		if (((unsigned int)handle.idx() < _vc.size()) && (handle.idx() >= 0))
		{
			return _vc[handle.idx()];
		} 
		else
		{
			return *((Vertex*)NULL);
		}
	}

	/** Retrieve Point stored in the TetraMesh.
	*	\param  handle the index of the Point to retrieve, must be valid. 
	*   \return Point which indicated by the handle
	*/
	Point& TetraMesh::handle_to_entity(PointHandle& handle)
	{
		if (((unsigned int)handle.idx() < _pc.size()) && handle.idx() >= 0)
		{
			return _pc[handle.idx()];
		} 
		else
		{
			return *((Point*)NULL);
		}
	}

	/** Retrieve Point stored in the TetraMesh.
	*	\param  handle the index of the Point to retrieve, must be valid. 
	*   \return Point which indicated by the handle
	*/
	const Point& TetraMesh::handle_to_entity(PointHandle& handle) const
	{
		if (((unsigned int)handle.idx() < _pc.size()) && handle.idx() >= 0)
		{
			return _pc[handle.idx()];
		} 
		else
		{
			return *((Point*)NULL);
		}
	}

	
	/** Retrieve Tetra stored in the TetraMesh, translate the iterator to the item
	 *	\param  iter the iterator of the Tetra to retrieve, its handle must be valid. 
	 *   \return Tetrahedron which indicated by the iterator
	 */
	Tetrahedron & TetraMesh::iter_to_entity(HedronIter & iter)
	{
		if ( ((unsigned int)iter.handle().idx() < _tetc.size()) && (iter.handle().idx() >= 0))
		{
			return _tetc[iter.handle().idx()];
		}
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve Vertex stored in the TetraMesh indicated by iterator
	*	\param  iter the iterator of the Vertex to retrieve, its handle must be valid. 
	*   \return Vertex which indicated by the iterator
	*/
	Vertex& TetraMesh::iter_to_entity(VertexIter& iter)
	{
		if (((unsigned int)iter.handle().idx() < _vc.size()) && (iter.handle().idx() >= 0))
		{
			return _vc[iter.handle().idx()];
		} 
		else
		{
			return *((Vertex*)NULL);
		}
	}

	/** Retrieve Point stored in the TetraMesh indicated by iterator
	*   \param iter the iterator of the point to retrieve, its handle must be valid
	*   \return Point which indicated by the iterator
	*/
	Point & TetraMesh::iter_to_entity(PointIter& iter)
	{
		if (((unsigned int)iter.handle().idx() < _vc.size()) && (iter.handle().idx() >= 0))
		{
			return _pc[iter.handle().idx()];
		} 
		else
		{
			return *((Point *)NULL);
		}

	}


	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	Tetrahedron& TetraMesh::next_tetrahedron(TetraHandle &handle)
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return handle_to_entity(next_valid_tetrahedon_handle(handle));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	TetraHandle TetraMesh::next_tetrahedron_handle(TetraHandle &handle)
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return next_valid_tetrahedon_handle(handle);
		else
			return TetraHandle(-1);
	}

	/** Retrieve next HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	HalfFaceTH& TetraMesh::next_half_face(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle_to_entity(handle).next_half_face_handle());
	}

	/** Retrieve next HalfFaceHandle stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	HalfFaceHandle TetraMesh::next_half_face_handle(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle).next_half_face_handle();
	}

	/** Retrieve next HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	HalfEdgeTH& TetraMesh::next_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).next_half_edge_handle());
	}

	/** Retrieve next HalfEdgeHandle stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	HalfEdgeHandle TetraMesh::next_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).next_half_edge_handle();
	}

	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	Tetrahedron& TetraMesh::next_tetrahedron(Tetrahedron & th_)
	{
		if ( (th_.handle().idx() <= _npolyhedron - 1) && (th_.handle().idx() >= 0))
			return handle_to_entity(next_valid_tetrahedon_handle(th_.handle()));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	TetraHandle TetraMesh::next_tetrahedron_handle(Tetrahedron & th_)
	{
		if ( (th_.handle().idx() <= _npolyhedron - 1) && (th_.handle().idx() >= 0))
			return next_valid_tetrahedon_handle(th_.handle());
		else
			return TetraHandle(-1);
	}

	/** Retrieve next HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	HalfFaceTH& TetraMesh::next_half_face(HalfFaceTH& hf_)
	{
		return handle_to_entity(hf_.next_half_face_handle());
	}

	/** Retrieve next HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	HalfFaceHandle TetraMesh::next_half_face_handle(HalfFaceTH& hf_)
	{
		return hf_.next_half_face_handle();
	}

	/** Retrieve next HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	HalfEdgeTH& TetraMesh::next_half_edge(HalfEdgeTH &he_)
	{
		return handle_to_entity(he_.next_half_edge_handle());
	}

	/** Retrieve next HalfEdgeHandle stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	HalfEdgeHandle TetraMesh::next_half_edge_handle(HalfEdgeTH &he_)
	{
		return he_.next_half_edge_handle();
	}

	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	const Tetrahedron& TetraMesh::next_tetrahedron(TetraHandle &handle) const
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return handle_to_entity(next_valid_tetrahedon_handle(handle));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve next HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	const HalfFaceTH & TetraMesh::next_half_face(HalfFaceHandle & handle) const
	{
		HalfFaceHandle hh; 
		hh = handle_to_entity(handle).next_half_face_handle();
		return handle_to_entity(hh);
	}

	/** Retrieve next HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	const HalfEdgeTH& TetraMesh::next_half_edge(HalfEdgeHandle & handle) const
	{
		HalfEdgeHandle hh;
		hh = handle_to_entity(handle).next_half_edge_handle();
		return handle_to_entity(hh);
	}

	/** Retrieve next Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its next, must be valid. 
	*   \return Tetrahedron which is the next of the handle.
	*/
	const Tetrahedron& TetraMesh::next_tetrahedron(Tetrahedron & th_) const
	{
		if ( (th_.handle().idx() < _npolyhedron - 1) && (th_.handle().idx() >= 0))
			return handle_to_entity(TetraHandle(th_.handle().idx() + 1));
		else if ( th_.handle().idx() == _npolyhedron - 1 )
			return handle_to_entity(TetraHandle(0));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve next HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its next, must be valid. 
	*   \return HalfFaceTH which is the next of the handle.
	*/
	const HalfFaceTH& TetraMesh::next_half_face(HalfFaceTH& hf_) const
	{
		return handle_to_entity(hf_.next_half_face_handle());
	}

	/** Retrieve next HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its next, must be valid. 
	*   \return HalfEdgeTH which is the next of the handle.
	*/
	const HalfEdgeTH& TetraMesh::next_half_edge(HalfEdgeTH & he_) const
	{
		return handle_to_entity(he_.next_half_edge_handle());
	}

	//-----------------------------------------------------------------------------------------------------------------


	/** Retrieve previous Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	Tetrahedron& TetraMesh::prev_tetrahedron(TetraHandle &handle)
	{
		if ( handle.idx() <= _npolyhedron - 1 && handle.idx() >= 0 )
			return handle_to_entity(prev_valid_tetrahedon_handle(handle));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve previous Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	TetraHandle TetraMesh::prev_tetrahedron_handle(TetraHandle &handle)
	{
		if ( handle.idx() <= _npolyhedron - 1 && handle.idx() >= 0 )
			return prev_valid_tetrahedon_handle(handle);
		else
			return TetraHandle(-1);
	}

	/** Retrieve previous HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return HalfFaceTH which is the previous of the handle.
	*/
	HalfFaceTH& TetraMesh::prev_half_face(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_face_handle());
	}

	/** Retrieve previous HalfFace  stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return HalfFaceTH which is the previous of the handle.
	*/
	HalfFaceHandle TetraMesh::prev_half_face_handle(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle).prev_half_face_handle();
	}

	/** Retrieve previous HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdgeTH which is the previous of the handle.
	*/
	HalfEdgeTH& TetraMesh::prev_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_edge_handle());
	}

	/** Retrieve previous HalfEdge  handle stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdgeTH which is the previous of the handle.
	*/
	HalfEdgeHandle TetraMesh::prev_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).prev_half_edge_handle();
	}

	/** Retrieve previous Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	Tetrahedron& TetraMesh::prev_tetrahedron(Tetrahedron &th_)
	{
		if ( (th_.handle().idx() <= _npolyhedron - 1) && (th_.handle().idx() >= 0) )
			return handle_to_entity(prev_valid_tetrahedon_handle(th_.handle()));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve previous Tetra handle stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	TetraHandle TetraMesh::prev_tetrahedron_handle(Tetrahedron &th_)
	{
		if ( (th_.handle().idx() <= _npolyhedron - 1) && (th_.handle().idx() >= 0) )
			return prev_valid_tetrahedon_handle(th_.handle());
		else
			return TetraHandle(-1);
	}

	/** Retrieve previous HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return HalfFaceTH which is the previous of the handle.
	*/
	HalfFaceTH& TetraMesh::prev_half_face(HalfFaceTH& hf_)
	{
		return handle_to_entity(hf_.prev_half_face_handle());
	}

	/** Retrieve previous HalfFace handle stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return halfFace handle which is the previous of the handle.
	*/
	HalfFaceHandle TetraMesh::prev_half_face_handle(HalfFaceTH& hf_)
	{
		return hf_.prev_half_face_handle();
	}


	/** Retrieve previous HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdgeTH which is the previous of the handle.
	*/
	HalfEdgeTH& TetraMesh::prev_half_edge(HalfEdgeTH &he_)
	{
		return handle_to_entity(he_.prev_half_edge_handle());
	}


	/** Retrieve previous HalfEdge handle stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdge handle which is the previous of the handle.
	*/
	HalfEdgeHandle TetraMesh::prev_half_edge_handle(HalfEdgeTH &he_)
	{
		return he_.prev_half_edge_handle();
	}



	/** Retrieve previous Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	const Tetrahedron& TetraMesh::prev_tetrahedron(TetraHandle &handle) const
	{
		if ( handle.idx() <= _npolyhedron - 1 && handle.idx() >= 0 )
			return handle_to_entity(prev_valid_tetrahedon_handle(handle));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve previous HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return HalfFaceTH which is the previous of the handle.
	*/
	const HalfFaceTH& TetraMesh::prev_half_face(HalfFaceHandle& handle) const
	{
		HalfFaceHandle hh;
		hh = handle_to_entity(handle).prev_half_face_handle(); 
		return handle_to_entity(hh);
	}

	/** Retrieve previous HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdgeTH which is the previous of the handle.
	*/
	const HalfEdgeTH& TetraMesh::prev_half_edge(HalfEdgeHandle &handle) const
	{
		HalfEdgeHandle hh;
		hh = handle_to_entity(handle).prev_half_edge_handle(); // handle_to_entity(handle).prev_half_face_handle(); 
		return handle_to_entity(hh);
	}

	/** Retrieve previous Tetra stored in the TetraMesh.
	*	\param  handle the index of the Tetra which we want to get its previous, must be valid. 
	*   \return Tetrahedron which is the previous of the handle.
	*/
	const Tetrahedron& TetraMesh::prev_tetrahedron(Tetrahedron &th_) const
	{
		if ( (th_.handle().idx() <= _npolyhedron - 1) && (th_.handle().idx() >= 0) )
			return handle_to_entity(prev_valid_tetrahedon_handle(th_.handle()));
		else
			return *((Tetrahedron*)NULL);
	}

	/** Retrieve previous HalfFaceTH stored in the TetraMesh.
	*	\param  handle the index of the HalfFaceTH which we want to get its previous, must be valid. 
	*   \return HalfFaceTH which is the previous of the handle.
	*/
	const HalfFaceTH& TetraMesh::prev_half_face(HalfFaceTH & hf_) const
	{
		return handle_to_entity(hf_.prev_half_face_handle());
	}


	/** Retrieve previous HalfEdgeTH stored in the TetraMesh.
	*	\param  handle the index of the HalfEdgeTH which we want to get its previous, must be valid. 
	*   \return HalfEdgeTH which is the previous of the handle.
	*/
	const HalfEdgeTH& TetraMesh::prev_half_edge(HalfEdgeTH & he_) const
	{
		return handle_to_entity(he_.prev_half_edge_handle());
	}


	//-----------------------------------------------------------------------------------------------------------------

	/** Access from HalfEdgeTH to its start Point handle.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return PointHandle which gives the start point of HalfEdgeTH.
	*/
	PointHandle TetraMesh::start_point_handle(HalfEdgeHandle &handle)
	{
		if ( ((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0))
		{			
			return _vc[_hec[handle.idx()].start_vertex_handle().idx()].point_handle();
		} 
		else
		{
			return PointHandle(-1);
		}
	}

	/** Access from HalfEdgeTH to its Point.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Point which is the start point of HalfEdgeTH.
	*/
	Point& TetraMesh::start_point(HalfEdgeHandle & he_)
	{
		return handle_to_entity(start_point_handle(he_));
	}

	/** Retrieve the vertex handle which the half edge is from
	 *  \param handle the handle of the half edge, must be valid
	 *  \return VertexHandle the vertex handle 
	 */
	VertexHandle TetraMesh::from_vertex_handle(HalfEdgeHandle hh_)
	{
		assert(hh_.is_valid());
		return handle_to_entity(hh_).start_vertex_handle();
	}

	/** Retrieve the vertex handle which the half edge is to
	*  \param handle the handle of the half edge, must be valid
	*  \return VertexHandle the vertex handle 
	*/
	VertexHandle TetraMesh::to_vertex_handle(HalfEdgeHandle hh_)
	{
		assert(hh_.is_valid());
		return next_half_edge(hh_).start_vertex_handle();
	}

	/** Access from point handle to its Point.
	*	\param  handle the index of the point , must be valid. 
	*   \return Point which is the start point of HalfEdgeTH.
	*/
	Point& TetraMesh::point(PointHandle & handle)
	{		
		return handle_to_entity(handle);
	}

	/** Access from vertex handle to its Point.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Point which is the start point of HalfEdgeTH.
	*/
	Point& TetraMesh::point(VertexHandle & handle)
	{
		return handle_to_entity(handle_to_entity(handle).point_handle());
	}


	/** Access from point iterator to its Point.
	*	\param  iter the iterator of the point, its handle must be valid. 
	*   \return Point which is the point of iterator .
	*/
	Point& TetraMesh::point(PointIter & iter)
	{		
		return handle_to_entity(iter.handle());
	}

	/** Access from vertex handle to its Point.
	*	\param  handle the index of the vertex, its handle must be valid. 
	*   \return Point which is the point of iterator.
	*/
	Point& TetraMesh::point(VertexIter & iter )
	{
		return handle_to_entity(handle_to_entity(iter.handle()).point_handle());
	}



	//-----------------------------------------------------------------------------------------------------------------
	/** Retrieve the first half face handle to from a tetrahedron handle
	*	\param  th_ the index of the tetrahedron, its handle must be valid. 
	*   \return HalfFaceHandle which is the handle of the half face.
	*/
	HalfFaceHandle TetraMesh::first_half_face_handle(TetraHandle th_)
	{
		return handle_to_entity(th_).first_half_face_handle();
	}

	/** Retrieve the first half edge handle to from a half handle
	*	\param  fh_ the index of the half face, its handle must be valid. 
	*   \return HalfEdgeHandle which is the handle of the half edge.
	*/
	HalfEdgeHandle TetraMesh::first_half_edge_handle(HalfFaceHandle fh_)
	{
		return handle_to_entity(fh_).first_half_edge_handle();
	}

	/** Retrieve the first half face handle to from a tetrahedron iterator
	*	\param  tIter the iterator of the tetrahedron, its handle must be valid. 
	*   \return HalfFaceHandle which is the handle of the half face.
	*/
	HalfFaceHandle TetraMesh::first_half_face_handle(HedronIter tIter_)
	{
		return iter_to_entity(tIter_).first_half_face_handle();
	}

	/** Retrieve the first vertex handle to from a half face handle
	*	\param  fh_ the handle of the half face, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle TetraMesh::first_vertex_handle(HalfFaceHandle fh_)
	{
		return handle_to_entity(handle_to_entity(fh_).first_half_edge_handle()).start_vertex_handle();
	}

	/** Retrieve the first vertex handle to from a hedron handle
	*	\param  hh_ the handle of the tetrahedron, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle TetraMesh::first_vertex_handle(PolyhedronHandle hh_)
	{
		return handle_to_entity(hh_).first_vertex_handle();
	}

	/** Retrieve the first vertex handle  from a hedron iter
	*	\param  tIter_ the iterator of the tetrahedron, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle TetraMesh::first_vertex_handle(HedronIter tIter_)
	{
		return handle_to_entity(tIter_.handle()).first_vertex_handle();
	}

	/** Retrieve the next vertex handle from a hedron handle 
	 *  \param hh_ the polyhedral handle
	 *  \param cvh_ the current vertex handle
	 *  \return VertexHandle which is the next vertex handle
	 */
	VertexHandle TetraMesh::next_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && hh_.is_valid());
		return VertexHandle ((cvh_.idx() + 1) % 4 + handle_to_entity(hh_).first_vertex_handle().idx()); //((int)(cvh_.idx() / 4)) * 4
	}

	/** Retrieve the next vertex handle from a hedron iterator 
	*  \param tIter_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle TetraMesh::next_vertex_handle(HedronIter tIter_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && tIter_.handle().is_valid());
		return VertexHandle((tIter_.handle().idx() + 1) % 4 + handle_to_entity(tIter_.handle()).first_vertex_handle().idx());//((int)(hh_.idx() / 4)) * 4
	}

	/** Retrieve the prev vertex handle from a hedron handle 
	*  \param hh_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle TetraMesh::prev_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && hh_.is_valid());
		return VertexHandle ((cvh_.idx() - 1) % 4 + handle_to_entity(hh_).first_vertex_handle().idx()); //((int)(cvh_.idx() / 4)) * 4
	}

	/** Retrieve the prev vertex handle from a hedron iterator 
	*  \param tIter_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle TetraMesh::prev_vertex_handle(HedronIter tIter_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && tIter_.handle().is_valid());
		return VertexHandle((tIter_.handle().idx() - 1) % 4 + handle_to_entity(tIter_.handle()).first_vertex_handle().idx());//((int)(hh_.idx() / 4)) * 4
	}

	//-----------------------------------------------------------------------------------------------------------------
	/** Retrieve the first half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::first_hedron_half_edge_handle(HedronIter & hIter_)
	{
		assert(hIter_.handle().is_valid());
		return HalfEdgeHandle(hIter_.handle().idx() * 12);
	}

	/** Retrieve the first half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::first_hedron_half_edge_handle(HedronHandle & hh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12);
	}

	extern const int TNEXT[4][4];
	/** Retrieve the first half face handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::first_vertex_half_face_handle(VertexHandle & vh_)
	{
		assert(vh_.is_valid());
		int idx;
		for (int i = 0; i < 4; i ++)
		{
			if (TNEXT[i][vh_%4] != -1)
			{
				idx = TNEXT[i][vh_%4];
				break;
			}
		}
		return HalfFaceHandle((vh_>>2)*4+idx);
	}
	
	/** Retrieve the first half face handle related to a vertex
	*   \param in vIter the vertex iterator
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::first_vertex_half_face_handle(VertexIter   & vIter_)
	{
		VertexHandle vh_ = vIter_.handle();
		assert(vh_.is_valid());
		int idx;
		for (int i = 0; i < 4; i ++)
		{
			if (TNEXT[i][vh_%4] != -1)
			{
				idx = TNEXT[i][vh_%4];
				break;
			}
		}
		return HalfFaceHandle((vh_>>2)*4+idx);
	}


	/** Retrieve the first half edge handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::first_vertex_half_edge_handle(VertexHandle & vh_)
	{
		assert(vh_.is_valid());
		return HalfEdgeHandle((vh_>>2)*12+vhecTH_[vh_%4]);
	}
		
	/** Retrieve the first half edge handle related to a vertex
	*   \param in vIter the vertex iterator
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::first_vertex_half_edge_handle(VertexIter   & vIter_)
	{
		assert(vIter_.handle().is_valid());
		return HalfEdgeHandle((vIter_.handle()>>2)*4+vhecTH_[vIter_.handle()%4]);
	}

	/** Retrieve the next half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::next_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_)
	{
		assert(hIter_.handle().is_valid() && heh_.is_valid());
		return HalfEdgeHandle(hIter_.handle() * 12 + (heh_.idx() + 1) % 12);
	}

	/** Retrieve the next half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::next_hedron_half_edge_handle(HedronHandle & hh_, HalfEdgeHandle & heh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12 + (heh_.idx() + 1) % 12);
	}

	/** Retrieve the previous half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::prev_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_)
	{
		assert(hIter_.handle().is_valid() && heh_.is_valid());
		return HalfEdgeHandle(hIter_.handle() * 12 + (heh_.idx() - 1) % 12);
	}

	/** Retrieve the previous half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle TetraMesh::prev_hedron_half_edge_handle(HedronHandle & hh_, HalfEdgeHandle & heh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12 + (heh_.idx() - 1) % 12);
	}

	/** Retrieve the next half face handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::next_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_)
	{
		assert(vh_.is_valid());
		return HalfFaceHandle((hfh_>>2)*4+TNEXT[hfh_%4][vh_%4]);
	}
	
	/** Retrieve the next half face handle related to a vertex
	*   \param in vIter_ the vertex iterator
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::next_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_)
	{
		assert(vIter_.handle().is_valid());
		return HalfFaceHandle((hfh_>>2)*4+TNEXT[hfh_%4][vIter_.handle()%4]);
	}

	/** Retrieve the prev half face handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::prev_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_)
	{
		assert(vh_.is_valid());
		return HalfFaceHandle((hfh_>>2)*4+TNEXT[TNEXT[hfh_%4][vh_%4]][vh_%4]);
	}

	/** Retrieve the prev half face handle related to a vertex
	*   \param in vIter_ the vertex iterator
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle TetraMesh::prev_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_)
	{
		assert(vIter_.handle().is_valid());
		return HalfFaceHandle((hfh_>>2)*4+TNEXT[TNEXT[hfh_%4][vIter_.handle()%4]][vIter_.handle()%4]);
	}
	

	//-----------------------------------------------------------------------------------------------------------------





	/** Get Radial HalfEdgeTH of a particular HalfEdgeTH.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Radial HalfEdgeTH.
	*/
	HalfEdgeTH& TetraMesh::radial_half_edge(HalfEdgeHandle &handle)
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0))
		{
			HalfEdgeTH he = handle_to_entity(handle);
			HalfFaceTH hf  = handle_to_entity(he.half_face_handle());
			if (hf.has_opposite_face())
			{
				PointHandle nextp = start_point_handle(he.next_half_edge_handle());
				HalfFaceTH opphf = opposite_half_face(hf.handle());
				HalfEdgeTH opphe = handle_to_entity(opphf.first_half_edge_handle());
				while( nextp.idx() != start_point_handle(opphe.handle()).idx() )
				{
					opphe = handle_to_entity(opphe.next_half_edge_handle());
				}
				return handle_to_entity(opphe.handle());
			}
			else
				return *((HalfEdgeTH*)NULL);
		} 
		else
		{
			return *((HalfEdgeTH*)NULL);
		}
	}

	/** Get Radial HalfEdge handle of a particular HalfEdgeTH.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Radial HalfEdgeTH.
	*/
	HalfEdgeHandle TetraMesh::radial_half_edge_handle(HalfEdgeHandle &handle)
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0))
		{
			HalfEdgeTH he = handle_to_entity(handle);
			HalfFaceTH hf  = handle_to_entity(he.half_face_handle());
			if (hf.has_opposite_face())
			{
				PointHandle nextp = start_point_handle(he.next_half_edge_handle());
				HalfFaceTH opphf = opposite_half_face(hf.handle());
				HalfEdgeTH opphe = handle_to_entity(opphf.first_half_edge_handle());
				while( nextp.idx() != start_point_handle(opphe.handle()).idx() )
				{
					opphe = handle_to_entity(opphe.next_half_edge_handle());
				}
				return opphe.handle();
			}
			else
				return HalfEdgeHandle(-1);
		} 
		else
		{
			return HalfEdgeHandle(-1);
		}
	}


	/** Check whether the HalfEdgeTH has its Radial HalfEdgeTH.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_radial_half_edge(HalfEdgeHandle &handle)
	{
		HalfEdgeTH he = handle_to_entity(handle);
		HalfFaceTH hf  = handle_to_entity(he.half_face_handle());
		return hf.has_opposite_face();
	}

	/** Get mate HalfEdgeTH of a particular HalfEdgeTH.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Mate HalfEdgeTH.
	*/
	HalfEdgeTH& TetraMesh::mate_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).mate_half_edge_handle());
	}

	/** Get mate HalfEdgeTH of a particular HalfEdgeTH.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Mate HalfEdgeTH.
	*/
	HalfEdgeHandle TetraMesh::mate_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).mate_half_edge_handle();
	}

	/** Get Opposite HalfFaceTH of a particular HalfFaceTH.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return Opposite HalfFaceTH.
	*/
	HalfFaceTH& TetraMesh::opposite_half_face(HalfFaceHandle& handle)
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0) )
		{
			HalfFaceTH hf = handle_to_entity(handle);
			if (hf.has_opposite_face())
				return handle_to_entity(hf.opposite_face_handle());
			else
				return *((HalfFaceTH*)NULL);
		}
		else
		{
			return *((HalfFaceTH*)NULL);
		}
	}

	/** Get Opposite HalfFace handle of a particular HalfFaceTH.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return Opposite HalfFaceTH.
	*/
	HalfFaceHandle TetraMesh::opposite_half_face_handle(HalfFaceHandle& handle)
	{
		if (((unsigned int)handle.idx() < _hec.size()) && (handle.idx() >= 0) )
		{
			HalfFaceTH hf = handle_to_entity(handle);
			if (hf.has_opposite_face())
				return hf.opposite_face_handle();
			else
				return HalfFaceHandle(-1);
		}
		else
		{
			return HalfFaceHandle(-1);
		}
	}

	//-----------------------------------------------------------------------------------------------------------------



	//-----------------------------------------------------------------------------------------------------------------

	/** Check whether the HalfFaceTH has its Opposite HalfFaceTH.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_opposite_half_face(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle).has_opposite_face();
	}

	/** Check whether the hedron handle has boundary half face.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_boundary_face(PolyhedronHandle & hh_)
	{
		HedronFaceIter hf_it;
		for (hf_it = hedron_face_iter(hh_); hf_it; ++ hf_it)
		{
			if (!has_opposite_half_face(hf_it.handle()))
			{
				return true;
			}
		}
		return false;
	}

	/** Check whether the HedronIter has boundary half face.
	*	\param  hIter_ iterator of tetrahedron, its handle must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_boundary_face(HedronIter & hIter_)
	{
		HedronFaceIter hf_it;
		for (hf_it = hedron_face_iter(hIter_); hf_it; ++ hf_it)
		{
			if (!has_opposite_half_face(hf_it.handle()))
			{
				return true;
			}
		}
		return false;
	}


	/** Get Tetrahedron of a particular Vertex.
	*	\param  handle the index of the Vertex, must be valid. 
	*   \return Tetrahedron.
	*/
	Tetrahedron& TetraMesh::get_tetrahedron(VertexHandle& handle)
	{
		return handle_to_entity(TetraHandle(handle.idx() >> 2));
	}

	/** Get Tetrahedron handle of a particular Vertex.
	*	\param  handle the index of the Vertex, must be valid. 
	*   \return Tetrahedron.
	*/
	TetraHandle TetraMesh::tetrahedron_handle(VertexHandle& handle)
	{
		return TetraHandle(handle.idx() >> 2);
	}

	/** Check whether the HalfFaceTH has particular Vertex.
	*	\param  hfHandle the index of the HalfFaceTH, must be valid. 
	*	\param  vHandle the index of the Vertex, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_vertex(HalfFaceHandle hfHandle, VertexHandle vHandle)
	{
		HalfFaceTH hf = handle_to_entity(hfHandle);
		HalfEdgeTH he = handle_to_entity(hf.first_half_edge_handle());
		for ( int i=0; i<3; i++ )
		{
			if ( he.start_vertex_handle() == vHandle )
				return true;
			he = next_half_edge(he.handle());
		}
		return false;
	}

	/** Check whether the HalfFaceTH has particular Point.
	*	\param  hfHandle the index of the HalfFaceTH, must be valid. 
	*	\param  pHandle the index of the Point, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::has_point(HalfFaceHandle hfHandle, PointHandle pHandle)
	{
		HalfFaceTH hf = handle_to_entity(hfHandle);
		HalfEdgeTH he = handle_to_entity(hf.first_half_edge_handle());
		for ( int i = 0; i<3; i++ )
		{
			if ( start_point_handle(he.handle()) == pHandle )
				return true;
			he = next_half_edge(he.handle());
		}
		return false;
	}

	/** Release all the topology stored in the TetraMesh class.
	*/
	void TetraMesh::release_mesh()
	{
		_hec.clear();
		_hfc.clear();
		_tetc.clear();
		_pc.clear();
		_vc.clear();
		_ophfc.clear();
		_npolyhedron = 0;
		_npoint = 0;
	}

	/** Add a Tetra to the TetrasMesh.
	*	\param  p1 point 1. 
	*	\param  p2 point 2.
	*	\param  p3 point 3.
	*	\param  p4 point 4.
	*   \return Tetrahandle(-1)  invalid parameter or already has this tetra.
	*			Tetrahandle(>=0) success.
	*/
	TetraHandle TetraMesh::add_tetrahedron(Point &p1, Point &p2, Point &p3, Point &p4)
	{
		using namespace std;
		int unused = -1;
		PointFindIF<Point> pif(p1);
		vector<int> pv;
		TetraHandle th;
		int i;
		int pidx;

		// check if the tetra is reverse
		Normal n[2];
		n[0] = (p3 - p2) % (p4 - p2);
		n[1] = p2 - p1;
		if ((n[0] | n[1]) < 0)
		{
			Point tp = p2;
			p2 = p4;
			p4 = tp;
		}

		// insert into empty location, just pop free_list and adjust pointHandle
		if ( _free_tetra_list.empty()==false )
		{
			th = pop_erase_from_free_tetra_list();
			i = th.idx();
			
			//reset halfface again to avoid invalid case
			for ( int j = 0; j < 4; j ++ )
			{
				_vc[i*4+j].set_next_vertex_handle(VertexHandle(-1));
				_hfc[i*4+j] = HalfFaceTH( i * 4 + j);
				_hec[(i*4+j)*3]   = HalfEdgeTH((i * 4 + j) * 3,	    i * 4 + (0 < j ? 0 : 1));
				_hec[(i*4+j)*3+1] = HalfEdgeTH((i * 4 + j) * 3 + 1, i * 4 + (1 < j ? 1 : 2));
				_hec[(i*4+j)*3+2] = HalfEdgeTH((i * 4 + j) * 3 + 2, i * 4 + (2 < j ? 2 : 3));
			}
		}
		// insert at tail
		else
		{
			i = _tetc.size();
			for ( int j = 0; j < 4; j ++ )
			{
				_hfc.push_back(HalfFaceTH( i * 4 + j));
				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3,	   i * 4 + (0 < j ? 0 : 1)));
				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3 + 1, i * 4 + (1 < j ? 1 : 2)));
				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3 + 2, i * 4 + (2 < j ? 2 : 3)));
				_vc.push_back(Vertex(-1));
			}
			_tetc.push_back(i);
			this->_npolyhedron++;
		}

		//find & add the vertex & point first
		// add new point
		if( (pidx=find_point(pif).idx())==-1 )
		{
			//add at tail
			if ( _free_point_list.empty()==true )
			{
				_pvc.push_back(VertexHandle(i*4));
				_pc.push_back(p1);
				_vc[i*4] = Vertex(i*4, _npoint);
				//_pvm[PointHandle(_npoint)].push_back(VertexHandle(i*4));
				this->_npoint++;
			}
			//add at free list location
			else
			{
				PointHandle ph = pop_erase_from_free_point_list();
				Point& p = handle_to_entity(ph);
				p = p1;
				_vc[i*4] = Vertex(i*4, ph.idx());
				_pvc[ph.idx()] = VertexHandle(i*4);
				//_pvm[ph].push_back(VertexHandle(i*4));
			}
			unused = 0;
		}
		// old point
		else
		{
			pv.push_back(pidx);
			_vc[i*4] = Vertex(i*4, pidx);
			_vc[i*4].set_next_vertex_handle(_pvc[pidx]);
			_pvc[pidx] = _vc[i*4].handle();
			//_pvm[PointHandle(pidx)].push_back(VertexHandle(i*4));
		}
		pif = p2;
		// add new point
		if((pidx=find_point(pif).idx())==-1)
		{
			//add at tail
			if ( _free_point_list.empty()==true )
			{
				_pvc.push_back(VertexHandle(i*4+1));
				_pc.push_back(p2);
				_vc[i*4+1] = Vertex(i*4+1, _npoint);
				//_pvm[PointHandle(_npoint)].push_back(VertexHandle(i*4+1));
				this->_npoint++;
			}
			//add at free list location
			else
			{
				PointHandle ph = pop_erase_from_free_point_list();
				Point& p = handle_to_entity(ph);
				p = p2;
				_vc[i*4+1] = Vertex(i*4+1, ph.idx());
				_pvc[ph.idx()] = VertexHandle(i*4+1);
				//_pvm[ph].push_back(VertexHandle(i*4+1));
			}
			unused = 1;
		}
		// old point
		else
		{
			pv.push_back(pidx);
			_vc[i*4+1] = Vertex(i*4+1, pidx);
			_vc[i*4+1].set_next_vertex_handle(_pvc[pidx]);
			_pvc[pidx] = _vc[i*4+1].handle();
			//_pvm[PointHandle(pidx)].push_back(VertexHandle(i*4+1));
		}
		pif = p3;
		// add new point
		if((pidx=find_point(pif).idx())==-1)
		{
			//add at tail
			if ( _free_point_list.empty()==true )
			{
				_pvc.push_back(VertexHandle(i*4+2));
				_pc.push_back(p3);
				_vc[i*4+2] = Vertex(i*4+2, _npoint);
				//_pvm[PointHandle(_npoint)].push_back(VertexHandle(i*4+2));
				this->_npoint++;
			}
			//add at free list location
			else
			{
				PointHandle ph = pop_erase_from_free_point_list();
				Point& p = handle_to_entity(ph);
				p = p3;
				_vc[i*4+2] = Vertex(i*4+2, ph.idx());
				_pvc[ph] = VertexHandle(i*4+2);
				//_pvm[ph].push_back(VertexHandle(i*4+2));
			}
			unused = 2;
		}
		// old point
		else
		{
			pv.push_back(pidx);
			_vc[i*4+2] = Vertex(i*4+2, pidx);
			_vc[i*4+2].set_next_vertex_handle(_pvc[pidx]);
			_pvc[pidx] = _vc[i*4+2].handle();
			//_pvm[PointHandle(pidx)].push_back(VertexHandle(i*4+2));
		}
		pif = p4;
		// add new point
		if((pidx=find_point(pif).idx())==-1)
		{
			//add at tail
			if ( _free_point_list.empty()==true )
			{
				_pvc.push_back(VertexHandle(i*4+3));
				_pc.push_back(p4);
				_vc[i*4+3] = Vertex(i*4+3, _npoint);
				//_pvm[PointHandle(_npoint)].push_back(VertexHandle(i*4+3));
				this->_npoint++;
			}
			//add at free list location
			else
			{
				PointHandle ph = pop_erase_from_free_point_list();
				Point& p = handle_to_entity(ph);
				p = p4;
				_vc[i*4+3] = Vertex(i*4+3, ph.idx());
				_pvc[ph] = VertexHandle(i*4+3);
				//_pvm[ph].push_back(VertexHandle(i*4+3));
			}
			unused = 3;
		}
		// old point
		else
		{
			pv.push_back(pidx);
			_vc[i*4+3] = Vertex(i*4+3, pidx);
			_vc[i*4+3].set_next_vertex_handle(_pvc[pidx]);
			_pvc[pidx] = _vc[i*4+3].handle();
			//_pvm[PointHandle(pidx)].push_back(VertexHandle(i*4+3));
		}

		//------------new add--------------//
		//adjust_tetrahedron(TetraHandle(i));  // have done before
		//------------new add--------------//

		int pointHandle[3];
		char buf[128];
		int v1,v2,v3;

		//there's not any half face share
		if ( pv.size() < 3 )
		{
			// set map relation 
			for ( int j = 0; j < 4; j++)
			{
				v1 = i * 4 + (0 < j ? 0 : 1);
				v2 = i * 4 + (1 < j ? 1 : 2);
				v3 = i * 4 + (2 < j ? 2 : 3);
				pointHandle[0] = _vc[v1].point_handle().idx();
				pointHandle[1] = _vc[v2].point_handle().idx();
				pointHandle[2] = _vc[v3].point_handle().idx();
				sort(pointHandle, pointHandle + 3);
				sprintf_s(buf,"%d %d %d", pointHandle[0], pointHandle[1], pointHandle[2]);
				_ophfc[buf] = i*4+j;
			}
			
			return TetraHandle(i);
		}
		// share 1 half face
		else if ( pv.size()==3 )
		{
			for ( int j=0; j<4; j++)
			{
				if ( j!=unused )
				{
					v1 = i*4 + (0<j ? 0 : 1);
					v2 = i*4 + (1<j ? 1 : 2);
					v3 = i*4 + (2<j ? 2 : 3);
					pointHandle[0] = _vc[v1].point_handle().idx();
					pointHandle[1] = _vc[v2].point_handle().idx();
					pointHandle[2] = _vc[v3].point_handle().idx();
					sort(pointHandle,pointHandle+3);
					sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
					_ophfc[buf] = i*4+j;
				} 
				else
				{
					sort(pv.begin(),pv.end());
					sprintf_s(buf,"%d %d %d",pv[0],pv[1],pv[2]);
					int idxOpHf = _ophfc[buf];
					if ( idxOpHf==-1 )
					{
						assert(false);
						cout<<"Error! There is more than one half face opposite to a particular half face."<<endl;
						
						if (i == _tetc.size()-1)
						{
							_pvc[_vc[i*4  ].point_handle().idx()] = _vc[i*4  ].next_vertex_handle();
							_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();
							_pvc[_vc[i*4+2].point_handle().idx()] = _vc[i*4+2].next_vertex_handle();
							_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();

							_tetc.pop_back();
							--_npolyhedron;
							_hfc.pop_back();_hfc.pop_back();
							_hfc.pop_back();_hfc.pop_back();
							for (int k=0;k<12;k++)
								_hec.pop_back();					
							_vc.pop_back();	_vc.pop_back();
							_vc.pop_back();	_vc.pop_back();
						}
						else
						{
							_pvc[_vc[i*4  ].point_handle().idx()] = _vc[i*4  ].next_vertex_handle();
							_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();
							_pvc[_vc[i*4+2].point_handle().idx()] = _vc[i*4+2].next_vertex_handle();
							_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();
							add_to_free_tetra_list(_tetc[i].handle());
						}

						for (int i = 0; i < 4-pv.size(); i++)
						{
							_pc.pop_back();
							_pvc.pop_back();
						}
						this->_npoint -= 4-pv.size();
						return TetraHandle(-1);
					}
					// reset map & opposite relation
					else
					{
						_hfc[idxOpHf].set_opposite_face_handle(HalfFaceHandle(i * 4 + j));
						_hfc[i * 4 + j].set_opposite_face_handle(HalfFaceHandle(idxOpHf));
						_ophfc[buf] = -1;
					}
				}
			}
			return TetraHandle(i);
		}
		//may share more than 1 half face
		else if ( pv.size() == 4 )
		{
			//dealing with half face shared with old tetra
			for ( int j=0; j<4; j++ )
			{
				v1 = i*4 + (0<j ? 0 : 1);	
				v2 = i*4 + (1<j ? 1 : 2);	
				v3 = i*4 + (2<j ? 2 : 3);
				pointHandle[0] = _vc[v1].point_handle().idx();
				pointHandle[1] = _vc[v2].point_handle().idx();
				pointHandle[2] = _vc[v3].point_handle().idx();
				sort(pointHandle,pointHandle+3);
				sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
				//find its opposite half face
				if ( _ophfc.find(buf)!=_ophfc.end() )
				{
					// reset map & opposite relation
					if ( _ophfc[buf]!=-1 )
					{
						int idxOpHf = _ophfc[buf];
						_hfc[idxOpHf].set_opposite_face_handle(HalfFaceHandle(i*4+j));
						_hfc[i*4+j].set_opposite_face_handle(HalfFaceHandle(idxOpHf));
						_ophfc[buf] = -1;
					} 
					else
					{
						assert(false);
						cout<<"Error! There is more than 1 half face opposite to another."<<endl;
						if (i == _tetc.size()-1)
						{
							_pvc[_vc[i*4  ].point_handle().idx()] = _vc[i*4  ].next_vertex_handle();
							_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();
							_pvc[_vc[i*4+2].point_handle().idx()] = _vc[i*4+2].next_vertex_handle();
							_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();

							_tetc.pop_back();
							--_npolyhedron;
							_hfc.pop_back();_hfc.pop_back();
							_hfc.pop_back();_hfc.pop_back();
							for (int k=0;k<12;k++)
								_hec.pop_back();					
							_vc.pop_back();	_vc.pop_back();
							_vc.pop_back();	_vc.pop_back();
						}
						else
						{
							_pvc[_vc[i*4  ].point_handle().idx()] = _vc[i*4  ].next_vertex_handle();
							_pvc[_vc[i*4+1].point_handle().idx()] = _vc[i*4+1].next_vertex_handle();
							_pvc[_vc[i*4+2].point_handle().idx()] = _vc[i*4+2].next_vertex_handle();
							_pvc[_vc[i*4+3].point_handle().idx()] = _vc[i*4+3].next_vertex_handle();
							add_to_free_tetra_list(_tetc[i].handle());
						}
						return TetraHandle(-1);
					}
				}
				else
				{
					_ophfc[buf] = i*4+j;
				}
			}
			return TetraHandle(i);
		}
		else
			return TetraHandle(-1);
	}

	
	/** Add a point to the TetrasMesh.
	*   \return new point's PointHandle
	*/
	PointHandle TetraMesh::add_point(Point p)
	{
		if (!_free_point_list.empty())
		{
			PointHandle ph = pop_erase_from_free_point_list();
			this->set_point(ph, p);
			_pvc[ph.idx()] = VertexHandle(-1);
			return ph;
		}

		_pc.push_back(p);
		_pvc.push_back(VertexHandle(-1));
		//_pvm[PointHandle(_npoint)].clear();
		++_npoint;
		return PointHandle(_npoint-1);
	}

	/** Erase a point of the TetraMesh.
	*	This function exchange the last point and the point to delete,
	*	then delete the last point and update the point and vertex map relation.
	*	\param  ph index of the point to erase. 
	*   \return -1 invalid parameter.
	*			 1 success.
	*/
	int TetraMesh::erase_point(PointHandle &ph)
	{
		//invalid handle
		if ( ph.idx() < 0 || ph.idx() >= _npoint || !is_valid(ph) )
			return -1;

		// add pointHandle to free_point_list.
		add_to_free_point_list(ph);

		//--------------------------new added---------------------------------------//
		//_pvc[ph.idx()] = VertexHandle(-1);
		//--------------------------end added---------------------------------------//

		return 1;
	}



	/** Erase a Tetra of the TetraMesh.
	*	This function invalid the tetra to erase by adding it to free_list,
	*	and update the opposite relation.
	*	\param  handle index of the tetra to erase. 
	*   \return next valid TetraHandle after the one erased,  success
	*			TetraHandle(-1)                            ,  no more tetra 
	*           TetraHandle(-1)                            ,  invalid input
	*/
	TetraHandle TetraMesh::erase_tetrahedron(TetraHandle &handle)
	{
		using namespace std;
		//invalid handle
		if ( handle.idx() < 0 || handle.idx() >= _npolyhedron )
		{
			std::cerr << "Error: TetraMesh::erase_tetrahedron/TetraHandle is out of size!" <<std::endl;
			return TetraHandle(-1);
		}

		//reset opposite HalfFaceTH relation & delete point if necessary
		Tetrahedron tetra = handle_to_entity(handle);
		HalfFaceTH hf = handle_to_entity(tetra.first_half_face_handle());
		set<PointHandle> ps;

		//reset opposite HalfFaceTH relation
		for ( int i=0; i<4; i++ )
		{
			char buf[128];
			PointHandle pointHandle[3];
			if ( hf.has_opposite_face() )
			{
				HalfFaceTH& ophf = handle_to_entity(hf.opposite_face_handle());
				HalfEdgeTH he = handle_to_entity(ophf.first_half_edge_handle());
				//reset opposite HalfFaceTH relation
				ophf.set_opposite_face_handle(HalfFaceHandle(-1));
				//reset opposite map
				pointHandle[0] = start_point_handle(he.handle());
				pointHandle[1] = start_point_handle(he.next_half_edge_handle());
				pointHandle[2] = start_point_handle(he.prev_half_edge_handle());
				sort(pointHandle, pointHandle+3);
				sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
				_ophfc[buf] = ophf.handle().idx();					
				ps.insert(pointHandle[0]);
				ps.insert(pointHandle[1]);
				ps.insert(pointHandle[2]);
			}
			else
			{
				HalfEdgeTH he= handle_to_entity(hf.first_half_edge_handle());
				//erase opposite map relation
				pointHandle[0] = start_point_handle(HalfEdgeHandle(he.handle()));
				pointHandle[1] = start_point_handle(HalfEdgeHandle(he.handle() + 1));
				pointHandle[2] = start_point_handle(HalfEdgeHandle(he.handle() + 2));
				sort(pointHandle, pointHandle + 3);
				sprintf_s(buf, "%d %d %d", pointHandle[0],pointHandle[1],pointHandle[2]);
				_ophfc.erase(buf);
			}
			hf = next_half_face(hf);
		}

		// update _pvc
		for (int i = 0; i < 4; i++)
		{
			Vertex v = handle_to_entity(VertexHandle(tetra.first_vertex_handle()+i));
			// if point related vertex is needed to be erasing, then update the _pvc[]
			if ( _pvc[v.point_handle().idx()] != v.handle())
				continue;

			VertexFindIF vif(v);
			VContainer::iterator loc = _vc.begin();;
			while( (loc=find_if(loc,_vc.begin()+tetra.first_vertex_handle(),vif))!=(_vc.begin()+tetra.first_vertex_handle()) )
			{
				if ( is_valid((*loc).handle()) )
					break;
				else 
				{
					loc++;
					if ( loc==(_vc.begin()+tetra.first_vertex_handle()) )
						break;
				}
			}
			if ( loc!=(_vc.begin()+tetra.first_vertex_handle()) )
			{
				// update _pvc
				if (_pvc[v.point_handle().idx()] == v.handle())
					_pvc[v.point_handle().idx()] = VertexHandle(loc-_vc.begin());
				continue;
			}

			loc = _vc.begin()+tetra.first_vertex_handle()+4;
			while( (loc=find_if(loc,_vc.end(),vif))!=_vc.end())
			{
				if ( is_valid((*loc).handle()) )
					break;
				else 
				{
					loc++;
					if ( loc==_vc.end() )
						break;
				}
			}
			if ( loc!=_vc.end() )
			{
				// update _pvc
				if (_pvc[v.point_handle().idx()] == v.handle())
					_pvc[v.point_handle().idx()] = VertexHandle(loc-_vc.begin());
				continue;
			}

			// the point needs to be deleted
			//if ( _pvc[v.point_handle().idx()] == v.handle())
		}

		//delete point
		if ( ps.size()<4 )
		{
			for ( int i=0; i<4; i++ )
			{
				Vertex v = handle_to_entity(VertexHandle(tetra.first_vertex_handle()+i));
				// This point might be delete. It doesn't belong to opposite tetra.
				if ( ps.find(v.point_handle())==ps.end() )
				{
					//check whether there are other vertex refers to the same point
					VertexFindIF vif(v);
					VContainer::iterator loc = _vc.begin();;
					while( (loc=find_if(loc,_vc.begin()+tetra.first_vertex_handle(),vif))!=(_vc.begin()+tetra.first_vertex_handle()) )
					{
						if ( is_valid((*loc).handle()) )
							break;
						else 
						{
							loc++;
							if ( loc==(_vc.begin()+tetra.first_vertex_handle()) )
								break;
						}
					}
					if ( loc!=(_vc.begin()+tetra.first_vertex_handle()) )
						continue;

					loc = _vc.begin()+tetra.first_vertex_handle()+4;
					while( (loc=find_if(loc,_vc.end(),vif))!=_vc.end())
					{
						if ( is_valid((*loc).handle()) )
							break;
						else 
						{
							loc++;
							if ( loc==_vc.end() )
								break;
						}
					}
					if ( loc!=_vc.end() )
						continue;

					// This point has no reference, should be delete. 
					// #1 add pointHandle to free_point_list.
					add_to_free_point_list(v.point_handle());

					//--------------------------new added---------------------------------------//
					_pvc[v.point_handle().idx()] = VertexHandle(-1);
					//--------------------------end added---------------------------------------//
				}
			}
		}

		//don't delete relation, but add them to free_list
		add_to_free_tetra_list(handle);

		return next_valid_tetrahedon_handle(handle);
	}

	/** Retrieve EdgeStar in TetraMesh.
	*	\param  eh handle of the HalfEdgeTH. 
	*	\param  tetraVec store the result of the EdgeStar. 
	*   \return number of the HalfEdgeTH in EdgeStar.
	*/
	int TetraMesh::edge_star(HalfEdgeHandle eh, std::vector<TetraHandle> &tetraVec)
	{
		if (eh.idx() < 0 || eh.idx() >= (int)_tetc.size() * 12)
		{
			std::cerr<< "Error: TetraMesh::edge_star/HalfEdgeHandle is out of size!" <<std::endl;
			return -1;
		}

		HalfEdgeTH beg = handle_to_entity(eh);
		HalfEdgeTH iter = beg;
		tetraVec.push_back(beg.hedron_handle());
		while ( has_radial_half_edge(iter.handle()) && radial_half_edge(iter.handle()).handle() != beg.mate_half_edge_handle() )
		{
			iter = radial_half_edge(iter.handle());
			iter = mate_half_edge(iter.handle());
			tetraVec.push_back(iter.hedron_handle());
		}

		if ( !has_radial_half_edge(iter.handle()) )
		{
			iter = mate_half_edge(beg.handle());
			while( has_radial_half_edge(iter.handle()) )
			{
				iter = radial_half_edge(iter.handle());
				iter = mate_half_edge(iter.handle());
				tetraVec.push_back(iter.hedron_handle());
			}
		}
		return tetraVec.size();
	}

	/** Retrieve VertexStar in TetraMesh.
	*	\param  vh handle of the Vertex. 
	*	\param  tetraVec store the result of the VertexStar. 
	*   \return number of the Vertex in VertexStar.
	*/
	int TetraMesh::vertex_star(VolumeMesh::VertexHandle vh, std::vector<TetraHandle> &tetraVec)
	{
		if (vh.idx() < 0 || vh.idx() >= (int)_vc.size())
		{
			std::cerr<< "Error: TetraMesh::vertex_star/VertexHandle is out of size!" <<std::endl;
			return -1;
		}

		Tetrahedron tetra = get_tetrahedron(vh);
		tetraVec.push_back(tetra.handle());
		for ( int i = 0; i<4; i++ )
		{
			if ( i==vh.idx()%4 )
				continue;

			HalfFaceTH hf = handle_to_entity(HalfFaceHandle(tetra.first_half_face_handle()+i));
			if ( hf.has_opposite_face() )
			{
				HalfFaceTH opphf = opposite_half_face(hf.handle());
				TetraHandle oppTetraHandle = opphf.hedron_handle();
				if(std::find(tetraVec.begin(),tetraVec.end(),oppTetraHandle)==tetraVec.end())
				{
					Tetrahedron oppTetra = handle_to_entity(oppTetraHandle);
					HedronVertexIter hv_it;
					for (hv_it = hedron_vertex_iter(oppTetra.handle()); hv_it; ++hv_it)
					{
						if (point_handle(hv_it.handle()) == point_handle(vh))
							break;
					}
					if (vertex_star(hv_it.handle(), tetraVec) < 0)
						return -1;
				}
			}
		}
		return tetraVec.size();
	}

	int TetraMesh::point_star(PointHandle ph, std::vector<TetraHandle> &tetraVec)
	{
		//if (ph.idx() < 0 || _pvm.find(ph) == _pvm.end() || !is_valid(ph))
		//	return 0;

		//std::vector<VertexHandle> v;
		//v = _pvm[ph];

		//for (unsigned int i = 0; i < v.size(); i++)
		//{
		//	if (is_valid(v[i]))
		//		tetraVec.push_back((get_tetrahedron(v[i])).handle());
		//}
		vertex_star(_pvc[ph.idx()], tetraVec);

		return tetraVec.size();
	}

	/** Retrieve Vertex half edge star in HexadMesh.
	*	\param  vh handle of the Vertex. 
	*	\param  halfEdgeVec store the result of the Vertex half edge Star. 
	*   \return number of the HalfEdge in halfEdgeVec.
	*/
	int TetraMesh::vertex_half_edge_star(VertexHandle vh, std::vector<HalfEdgeHandle> &halfEdgeVec)
	{
		std::vector<TetraHandle> tetraVec;
		VertexFindIF vif(handle_to_entity(vh));
		std::vector<Vertex>::iterator iter;
		VertexHalfEdgeIter vhe_it;
		Vertex currV(vh);
		vertex_star(vh, tetraVec);
		for (unsigned int i = 0; i < tetraVec.size(); i++)
		{
			if (!is_valid(tetraVec[i]))
				continue;
			iter = std::find_if(_vc.begin()+tetraVec[i]*4, _vc.begin()+tetraVec[i]*4+4,vif);
			currV = *iter;
			for (vhe_it = vertex_half_edge_iter(currV.handle()); vhe_it; ++ vhe_it)
				halfEdgeVec.push_back(vhe_it.handle());
		}
		return halfEdgeVec.size();
	}

	/** Retrieve Vertex half face star in HexadMesh.
	*	\param  vh handle of the Vertex. 
	*	\param  halfFaceVec store the result of the Vertex half face Star. 
	*   \return number of the HalfFace in halfFaceVec.
	*/
	int TetraMesh::vertex_half_face_star(VertexHandle vh, std::vector<HalfFaceHandle> &halfFaceVec)
	{
		std::vector<TetraHandle> tetraVec;
		VertexFindIF vif(handle_to_entity(vh));
		std::vector<Vertex>::iterator iter;
		VertexHalfFaceIter vhf_it;
		Vertex currV(vh);
		vertex_star(vh, tetraVec);
		for (unsigned int i = 0; i < tetraVec.size(); i++)
		{
			if (!is_valid(tetraVec[i]))
				continue;
			iter = std::find_if(_vc.begin()+tetraVec[i]*4, _vc.begin()+tetraVec[i]*4+4,vif);
			//assert(iter == NULL || iter == _vc.begin()+tetraVec[i]*4+4);
			currV = *iter;
			for (vhf_it = vertex_half_face_iter(currV.handle()); vhf_it; ++ vhf_it)
				halfFaceVec.push_back(vhf_it.handle());
		}
		return halfFaceVec.size();
	}

	
	/** Get half edge handle with a particular face and from vertex.
	*	\param  fh the index of the half face, must be valid. 
	*	\param  vh the index of the vertex, must be valid. 
	*   \return half edge handle.
	*/
	HalfEdgeHandle & TetraMesh::half_edge_handle(HalfFaceHandle & fh, VertexHandle & vh)
	{
		HalfEdgeTH he = handle_to_entity(handle_to_entity(fh).first_half_edge_handle());
		for ( int i = 0; i < 3; i ++ )
		{
			if ( he.start_vertex_handle() == vh )
				return he.handle();
			he = next_half_edge(he.handle());
		}
		return HalfEdgeHandle(-1);
	}

	/** Get half face handle with a particular half edge handle.
	*	\param  handle the index of the half edge, must be valid. 
	*   \return half face handle.
	*/
	HalfFaceHandle & TetraMesh::half_face_handle(HalfEdgeHandle & handle)
	{
		return handle_to_entity(handle).half_face_handle();
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** Build the boundary vertices of the mesh
	*/
	void TetraMesh::update_boundary_point()
	{
		HFContainer::iterator it;
		HalfFaceVertexIter hf_it;
		_brv.resize(_pc.size(), false);
		for (it = _hfc.begin(); it != _hfc.end(); ++ it)
		{
			if (is_boundary((*it).handle()))
			{
				for (hf_it = half_face_vertex_iter((*it).handle()); hf_it; ++ hf_it)
				{
					//_bv.insert(handle_to_entity(hf_it.handle()).point_handle());
					_brv[handle_to_entity(hf_it.handle()).point_handle().idx()] = true;
				}
			}
		}
	}

	/** build the opposite half face
	*   \return bool . true for success, false for not.
	*/
	bool TetraMesh::build_opposite_half_face()
	{
		int halfFaceNum = _npolyhedron * 4;
		int pointHandle[3];
		char buf[256];
		for ( int i = 0; i < halfFaceNum; i++ )
		{
			//Get 3 vertex of the halfFace
			HalfFaceTH &hf = _hfc[i];
			int v1, v2, v3;
			v1 = (i & ~3) + (0 < (i & 3) ? 0 : 1);	//i&~3 == (i>>2)<<2 == (i/4)*4
			v2 = (i & ~3) + (1 < (i & 3) ? 1 : 2);	//i%4 == i&3
			v3 = (i & ~3) + (2 < (i & 3) ? 2 : 3);

			pointHandle[0] = _vc[v1].point_handle().idx();
			pointHandle[1] = _vc[v2].point_handle().idx();
			pointHandle[2] = _vc[v3].point_handle().idx();
			//sort 3 vertex
			std::sort(pointHandle, pointHandle + 3);
			sprintf_s(buf, "%d %d %d", pointHandle[0], pointHandle[1], pointHandle[2]);
			if ( _ophfc.find(buf) != _ophfc.end())
			{
				if (_ophfc[buf] != -1)
				{
					int idxOpHf = _ophfc[buf];
					hf.set_opposite_face_handle(HalfFaceHandle(idxOpHf));
					_hfc[idxOpHf].set_opposite_face_handle(HalfFaceHandle(i));
					_ophfc[buf] = -1;
				}
				else
				{
					std::cerr << "Error! There is more than 1 half face opposite to a particular halfface." << std::endl;
					return false;
				}
			} 
			else
			{
				_ophfc[buf] = i;
			}
		}	
		update_boundary_point();	
		return true;
	}

	/** Build the topology of the mesh
	*  \return bool . true for success, false for not
	*/
	bool TetraMesh::build_topology()
	{
		int hfIdx, heIdx;
		_npoint = _pc.size();
		_npolyhedron = _vc.size() / 4;
		
		for ( int i = 0; i < _npolyhedron; i++ )
		{
			hfIdx = i * 4;
			heIdx = i * 12;
			_tetc.push_back(Tetrahedron(i));
			Normal n[2];
			n[0] = (_pc[_vc[i * 4 + 2].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]) % 
				   (_pc[_vc[i * 4 + 3].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]);
			n[1] = _pc[_vc[i * 4 + 1].point_handle().idx()] - _pc[_vc[i * 4].point_handle().idx()];

			if ((n[0] | n[1]) < 0)
			{
				std::swap(_vc[i * 4 + 1], _vc[i * 4 + 3]);
			}		

			for ( int j = 0; j < 4; j ++ )
			{
				_hfc.push_back(HalfFaceTH(i * 4 + j));

				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3,	   i * 4 + (0 < j ? 0 : 1)));
				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3 + 1, i * 4 + (1 < j ? 1 : 2)));
				_hec.push_back(HalfEdgeTH((i * 4 + j) * 3 + 2, i * 4 + (2 < j ? 2 : 3)));
			}
		}

		/**
		* build the connection between the points and the vertices & the mapping from points to vertices
		*/		
		_pvc.resize(_npoint, VertexHandle(-1));
		for (unsigned int i = 0; i < _vc.size(); i ++)
		{
			if (_pvc[_vc[i].point_handle().idx()].idx() != -1)
				_vc[i].set_next_vertex_handle(_pvc[_vc[i].point_handle().idx()]);
			_pvc[_vc[i].point_handle().idx()] = _vc[i].handle();
			//_pvm[_vc[i].point_handle()].push_back(_vc[i].handle());
		}
		
		return build_opposite_half_face();
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** update the face normals
	*/
	void TetraMesh::update_face_normals()
	{
		assert(_rfn);

		_fn.resize(_hfc.size(), Vec3d(1, 0, 0));
		HFContainer::iterator it;
		HalfFaceVertexIter hfv_iter;
		FaceHalfedgeIter fhe_it;
		Normal v[3];
		Normal n;
		//unsigned int i;
		for (it = _hfc.begin(); it != _hfc.end(); ++ it)
		{
			fhe_it = face_half_edge_iter((*it).handle());
			v[0] = start_point(fhe_it.handle());
			++ fhe_it;
			v[1] = start_point(fhe_it.handle());
			++ fhe_it;
			v[2] = start_point(fhe_it.handle());
			//i = 0;
			//for (hfv_iter = half_face_vertex_iter((*it).handle()); hfv_iter; ++ hfv_iter)
			//{
			//	v[i] = point(hfv_iter.handle());
			//	++ i;
			//}

			n = (v[1] - v[0]) % (v[2] - v[0]);
			n.normalize();

			_fn[(*it).handle().idx()] = n;
		}
	}

	/** update the face normal
	*/
	Vec3d TetraMesh::update_face_normal(HalfFaceHandle & hfh_)
	{
		//assert(_rfn);
		if (_fn.empty())
		{
			_fn.resize(_hfc.size(), Vec3d(1, 0, 0));
		}
		HalfFaceVertexIter hfv_iter;
		Vec3d v[3];
		Vec3d n;
		unsigned int i;
			i = 0;
			for (hfv_iter = half_face_vertex_iter(hfh_); hfv_iter; ++ hfv_iter)
			{
				v[i] = point(hfv_iter.handle());
				++ i;
			}

			n = (v[1] - v[0]) % (v[2] - v[0]);
			n.normalize();

			_fn[hfh_.idx()] = n;
			return _fn[hfh_.idx()];
	}

	/** return the face normal
	*/
	TetraMesh::Normal TetraMesh::normal(HalfFaceHandle &  hfh_)
	{
		return _fn[hfh_.idx()];
	}

	/** return the face normal
	*/
	TetraMesh::Normal TetraMesh::normal(const HalfFaceHandle &  hfh_) const
	{
		return _fn[hfh_.idx()];
	}

	/** Check whether the HalfFaceTH is boundary half face.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(HalfFaceHandle & hfh_)
	{
		return !handle_to_entity(hfh_).has_opposite_face();
	}

	/** Check whether the vertex handle  is boundary vertex.
	*	\param  vh_ the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(VertexHandle & vh_)
	{
		return _brv[handle_to_entity(vh_).point_handle().idx()];
		//return !(_bv.find(handle_to_entity(vh_).point_handle()) == _bv.end());
	}

	/** Check whether the point handle  is boundary point.
	*	\param  vh_ the index of the PointHandle, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(PointHandle & ph_)
	{
		return _brv[ph_.idx()];
	}

	/** Check whether the hedron is boundary hedron, if one of its vertices is boundary point, it is a boundary hedron
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(PolyhedronHandle & hh_)
	{
		HedronVertexIter hv_it;
		for (hv_it = hedron_vertex_iter(hh_); hv_it; ++ hv_it)
		{
			if (is_boundary(hv_it.handle()))
			{
				return true;
			}
		}
		return false;
	}

	/** Check whether the hedron is boundary hedron
	*	\param  handle the index of the HalfFaceTH, must be valid.
	*   \param type_ the boundary type. 0x0001 for half face boundary while 0x0002 for only vertex boundary 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(HedronHandle & hh_, unsigned int type_)
	{
		HedronFaceIter hf_it;
		bool isBoundary;
		isBoundary = false;
		for (hf_it = hedron_face_iter(hh_); hf_it; ++ hf_it)
		{
			if (is_boundary(hf_it.handle()))
			{
				isBoundary = true;
				break;
			}
		}
		if (type_ & 0x0001)
		{
			return isBoundary;
		}
		else if (type_ & 0x0002)
		{
			if (isBoundary)
			{
				/*
				* is face boundary, not only vertex boundary
				*/
				return false;
			}
			else
			{
				return is_boundary(hh_);
			}
		}

		return false;
	}

	/** Check whether the half edge is boundary half edge, if one of its star half face is boundary, it is a boundary half edge
	*	\param  handle the index of the HalfEdgeHandle, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_boundary(HalfEdgeHandle &heh_)
	{
		Point p_to, p_from, p1, p2;
		std::vector<HedronHandle> edgeStar;
		HedronFaceIter hf_it;
		FaceHalfedgeIter fe_it;
		p_from = point(from_vertex_handle(heh_));
		p_to = point(to_vertex_handle(heh_));

		edge_star(heh_, edgeStar);

		for (unsigned int i = 0; i < edgeStar.size(); i++)
		{
			for (hf_it = hedron_face_iter(edgeStar[i]); hf_it; ++ hf_it)
			{
				for (fe_it = face_half_edge_iter(hf_it.handle()); fe_it; ++fe_it)
				{
					p1 = point(from_vertex_handle(fe_it.handle()));
					p2 = point(to_vertex_handle(fe_it.handle()));
					if (((p1 == p_from && p2 == p_to) || (p1 == p_to && p2 == p_from)) && is_boundary(hf_it.handle()))
					{
						return true;
					}
				}
			}
		}
		return false;
	}

	/** Check whether the tetra is reverse
	*	\param  handle the index of tetra, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool TetraMesh::is_reverse(TetraHandle & th)
	{
		if (th.idx() < 0 || th.idx() >= (int)_tetc.size())
		{
			std::cerr << "Error: TetraMesh::is_reverse/TetraHandle is out of size!" << std::endl;
			return true;
		}

		int i = th.idx();
		Normal n[2];
		n[0] = (_pc[_vc[i * 4 + 2].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]) % 
			(_pc[_vc[i * 4 + 3].point_handle().idx()] - _pc[_vc[i * 4 + 1].point_handle().idx()]);
		n[1] = _pc[_vc[i * 4 + 1].point_handle().idx()] - _pc[_vc[i * 4].point_handle().idx()];

		if ((n[0] | n[1]) < 0)
		{
			return true;
		}
		return false;
	}

	/** contract an edge of the TetraMesh.
	*	This function contract the edge to a new point, the midpoint of the edge's two endpoints
	*	then adjust the topology the TetraMesh, at last delete the edge and Tetras shared the edge.
	*	\param heh index of the half edge to contract. 
	*   \param ph_new index of the new point
	*   \return -1 invalid parameter.
	*			 1 success.
	*/
	int TetraMesh::edge_contract(HalfEdgeHandle &he_, PointHandle &ph_new)
	{
		//invalid handle
		if ( he_.idx() < 0 || he_.idx() >= (int)_hec.size() )
			return -1;

		int epn;
		VertexHandle vf, vt;
		Point midpoint, pf_, pt_, pf, pt;
		PointHandle pht, phf; 
		std::vector<HedronHandle> edgeStarHH;
		std::vector<TetraHandle> vertexStar;
		HedronHalfEdgeIter hhe_it;
		// get midpoint and endpoints of the edge
		phf = point_handle(from_vertex_handle(he_));
		pht = point_handle(to_vertex_handle(he_));
		pf_ = point(from_vertex_handle(he_));
		pt_ = point(to_vertex_handle(he_));
		midpoint = (pf_ + pt_)/ 2.0;
		// get tetras which containing the edge
		edge_star(he_, edgeStarHH);
		// get tetras which containing the to_point
		vertex_star(to_vertex_handle(he_), vertexStar);

		// recovering data saving
		rec_ec_from_point = pf_;
		rec_ec_to_point = pt_;
		rec_ec_from_point_handle = phf;
		rec_ec_to_point_handle = pht;
		rec_ec_to_point_vertex_handle = _pvc[pht.idx()];
		rec_ec_erase_point.clear();
		rec_ec_erase_tetra.clear();
		rec_ec_erase_tetra = edgeStarHH;
		rec_ec_to_point_vertex_container.clear();

		// adjusting topology relation in the TetraMesh
		std::vector<VertexHandle>::iterator vec_it;
		for (unsigned int i = 0; i < edgeStarHH.size(); i++)
		{
			// find the contracting halfedge in the tetra & get endpoints
			for (hhe_it = hedron_half_edge_iter(edgeStarHH[i]); hhe_it; ++hhe_it)
			{
				vf = from_vertex_handle(hhe_it.handle());
				vt = to_vertex_handle(hhe_it.handle());
				pf = point(vf);
				pt = point(vt);
				if (pf == pf_ && pt == pt_)
					break;
			}
			if (pf != pf_ || pt != pt_)
				return -1;

			// get halfFaces opposite to the endpoints of the contracting edge
			HalfFaceTH hff = _hfc[from_vertex_handle(hhe_it.handle())];
			HalfFaceTH hft = _hfc[to_vertex_handle(hhe_it.handle())];
			// get opposite halfFaces of hff and hft
			HalfFaceHandle ophff = hff.opposite_face_handle();
			HalfFaceHandle ophft = hft.opposite_face_handle();
			
			char buf1[128], buf2[128];
			PointHandle pointHandle[3];
			HalfEdgeTH he_e = handle_to_entity(hff.first_half_edge_handle());
			pointHandle[0] = start_point_handle(he_e.handle());
			pointHandle[1] = start_point_handle(he_e.next_half_edge_handle());
			pointHandle[2] = start_point_handle(he_e.prev_half_edge_handle());
			std::sort(pointHandle, pointHandle+3);
			sprintf_s(buf1,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);

			he_e = handle_to_entity(hft.first_half_edge_handle());
			pointHandle[0] = start_point_handle(he_e.handle());
			pointHandle[1] = start_point_handle(he_e.next_half_edge_handle());
			pointHandle[2] = start_point_handle(he_e.prev_half_edge_handle());
			std::sort(pointHandle, pointHandle+3);
			sprintf_s(buf2,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);

			// erase tetra
			epn = _free_point_list.size();
			erase_tetrahedron(edgeStarHH[i]);
			// ¼ÇÂ¼É¾³ýµÄµã rec_ec_erase_point_handle
			for (int i = 0; i < (int)_free_point_list.size() - epn; i ++)
				rec_ec_erase_point.push_back(_free_point_list[(int)_free_point_list.size()-1-i]);

			// modify topology(reset opposite halfFace)
			if (ophff.idx() != -1)
				handle_to_entity(ophff).set_opposite_face_handle(ophft);
			if (ophft.idx() != -1)
				handle_to_entity(ophft).set_opposite_face_handle(ophff);

			// update opposite halfFace flag set _opphf[]
			_ophfc.erase(buf1);
			_ophfc[buf2] = -1;
		}

		//----------------------edge contracting-------------------------//
		// move from_point to the midpoint
		ph_new = phf;
		_pc[ph_new.idx()] = midpoint;

		HedronVertexIter hv_it;
		for (unsigned int i = 0; i < vertexStar.size(); i++)
		{
			if (!is_valid(vertexStar[i]))
				continue;
			for (hv_it = hedron_vertex_iter(vertexStar[i]); hv_it; ++hv_it)
			{
				if (point_handle(hv_it.handle()) == pht)
				{
					rec_ec_to_point_vertex_container.push_back(hv_it.handle());
					handle_to_entity(hv_it.handle()).set_point_handle(phf);
				}
			}
		}

		// move to_point to the midpoint(delete to_point & set vertex's PointHandle to be midpoint's handle)
		//std::vector<VertexHandle> mvh = _pvm[point_handle(to_vertex_handle(he_))];

		// delete to_point
		if (erase_point(pht) < 0)
			return -1;
		rec_ec_erase_point.push_back(pht);

		// add its vertex mapping to midpoint
		//for (unsigned int i = 0; i < mvh.size(); i++)
		//{
			//_pvm[ph_new].push_back(mvh[i]);
			//handle_to_entity(mvh[i]).set_point_handle(ph_new);
		//}
		//------------------- end edge contracting------------------------//
		return 1;
	}

	// clean _free_tetra_list and _free_point_list & filling the vacant
	void TetraMesh::clean_garbage()
	{
		// clean _free_tetra_list
		for (unsigned int i = 0; i < _free_tetra_list.size(); i++)
		{
			if (_free_tetra_list[i].idx() != _npolyhedron - 1)
			{
				Tetrahedron lastTetra = handle_to_entity(TetraHandle(_npolyhedron - 1));
				HalfFaceTH lastHf = handle_to_entity(lastTetra.first_half_face_handle());
				HalfFaceTH hf = handle_to_entity(handle_to_entity(_free_tetra_list[i]).first_half_face_handle());
				for ( int j=0; j<4; j++ )
				{
					//reset last HalfFaceTH's opposite relation to new position
					if ( lastHf.has_opposite_face() )
					{
						HalfFaceTH& ophf = handle_to_entity(lastHf.opposite_face_handle());
						ophf.set_opposite_face_handle(HalfFaceHandle(_free_tetra_list[i].idx()*4+j));
						_hfc[hf.handle()].set_opposite_face_handle(ophf.handle());
					}
					else
					{
						_hfc[hf.handle()].set_opposite_face_handle(HalfFaceHandle(-1));

						// reset _ophfc' Map relation for the last Tetra
						char buf[128];
						PointHandle pointHandle[3];
						HalfEdgeHandle lastHe = lastHf.first_half_edge_handle();
						pointHandle[0] = start_point_handle(HalfEdgeHandle(lastHe));
						pointHandle[1] = start_point_handle(HalfEdgeHandle(lastHe + 1));
						pointHandle[2] = start_point_handle(HalfEdgeHandle(lastHe + 2));
						std::sort(pointHandle, pointHandle+3);
						sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
						_ophfc[buf] = _free_tetra_list[i]*4+j;
					}
					hf = next_half_face(hf);
					lastHf = next_half_face(lastHf);

					//reset point handle 
					Vertex lastV = handle_to_entity(VertexHandle(lastTetra.first_vertex_handle().idx()+j));
					Vertex &v = handle_to_entity(VertexHandle((handle_to_entity(_free_tetra_list[i]).first_vertex_handle()).idx()+j));
					v.set_point_handle(lastV.point_handle());
				}
			}
			//delete all the topology relation
			_tetc.pop_back();
			--_npolyhedron;
			_hfc.pop_back();_hfc.pop_back();
			_hfc.pop_back();_hfc.pop_back();
			_vc.pop_back();_vc.pop_back();
			_vc.pop_back();_vc.pop_back();
			for ( int j=0;j<12;j++ )
				_hec.pop_back();
		}

		// clean _free_point_list
		for (unsigned int i = 0; i < _free_point_list.size(); i++)
		{
			if (_free_point_list[i].idx() != _npoint - 1)
			{
				for (unsigned int j = 0; j < _hfc.size(); j++ )
				{
					if ( has_point(HalfFaceHandle(j),PointHandle(_npoint-1)) )
					{
						char buf[128];
						int tmpHf;

						//Revise opposite mapping relation
						PointHandle pointHandle[3];
						pointHandle[0] = start_point_handle(HalfEdgeHandle(j*3));
						pointHandle[1] = start_point_handle(HalfEdgeHandle(j*3+1));
						pointHandle[2] = start_point_handle(HalfEdgeHandle(j*3+2));
						std::sort(pointHandle, pointHandle+3);
						sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
						if ( _ophfc.find(buf)!=_ophfc.end() )
						{
							tmpHf = _ophfc[buf];
							_ophfc.erase(buf);

							pointHandle[2] = _free_point_list[i];
							std::sort(pointHandle, pointHandle+3);
							sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
							_ophfc[buf] = tmpHf;
						}
					}
				}
				//Revise PointHandle in Vertex
				VertexFindIF vif(Vertex(_npoint-1));
				std::vector<Vertex>::iterator iter;
				while ( (iter=find_if(_vc.begin(),_vc.end(),vif))!=_vc.end() )
					(*iter).set_point_handle(_free_point_list[i]);
				//--------------------------new added---------------------------------------//
				_pvc[_free_point_list[i].idx()] = _pvc[_npoint - 1];
				_pvm[_free_point_list[i]] = _pvm[PointHandle(_npoint - 1)];
				_pvc.pop_back();
				_pvm.erase(_free_point_list[i]);
				//--------------------------end added---------------------------------------//
			}
			//------------------------delete point------------------------//
			_pc[_free_point_list[i]] = _pc[_npoint - 1];
			_pc.pop_back();
			_npoint--;
			//-------------------------end delete-------------------------//
		}
	}
	
	PointHandle TetraMesh::find_point(PointFindIF<Point>& pif)
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


	/** operation recover
	*   return : 1 succeed  
	*            0 failed
	*/
	int TetraMesh::recover_edge_contract()
	{
		// check if the memory space is still available
		std::vector<TetraHandle>::iterator _free_tetra_list_iter;
		std::vector<PointHandle>::iterator _free_point_list_iter;
		for (unsigned int i = 0; i < rec_ec_erase_point.size(); i++)
			if ((_free_point_list_iter = find(_free_point_list.begin(), _free_point_list.end(), rec_ec_erase_point[i])) == _free_point_list.end())
				return 0;

		for (unsigned int i = 0; i < rec_ec_erase_tetra.size(); i++)
			if ((_free_tetra_list_iter = find(_free_tetra_list.begin(), _free_tetra_list.end(), rec_ec_erase_tetra[i])) == _free_tetra_list.end())
				return 0;

		// recover point
		//_free_point_list.erase(_free_point_list_iter);
		for (unsigned int i = 0; i < rec_ec_erase_point.size(); i++)
			if ((_free_point_list_iter = find(_free_point_list.begin(), _free_point_list.end(), rec_ec_erase_point[i])) != _free_point_list.end())
				_free_point_list.erase(_free_point_list_iter);
		// recover tetrahedrons
		for (unsigned int i = 0; i < rec_ec_erase_tetra.size(); i++)
			if ((_free_tetra_list_iter = find(_free_tetra_list.begin(), _free_tetra_list.end(), rec_ec_erase_tetra[i])) != _free_tetra_list.end())
				_free_tetra_list.erase(_free_tetra_list_iter);

		// recover from_point and to_point's coordinate
		_pc[rec_ec_from_point_handle.idx()] = rec_ec_from_point;
		_pc[rec_ec_to_point_handle.idx()] = rec_ec_to_point;

		// recover to_point vertex map
		if (rec_ec_to_point_vertex_container.size())
		{
			_pvc[rec_ec_to_point_handle.idx()] = rec_ec_to_point_vertex_container[0];
			for (unsigned int i = 0; i < rec_ec_to_point_vertex_container.size(); i++)
			{
				handle_to_entity(rec_ec_to_point_vertex_container[i]).set_point_handle(rec_ec_to_point_handle);
			}
		}
		else
			_pvc[rec_ec_to_point_handle.idx()] = rec_ec_to_point_vertex_handle;

		// recover tetrahedrons
		HedronHalfEdgeIter hhe_it;
		char buf[128];
		PointHandle pointHandle[3];
		for (unsigned int i = 0; i < rec_ec_erase_tetra.size(); i++)
		{
			// find the contracting halfedge in the tetra & get endpoints
			for (hhe_it = hedron_half_edge_iter(rec_ec_erase_tetra[i]); hhe_it; ++hhe_it)
			{
				if (point(from_vertex_handle(hhe_it.handle())) == rec_ec_from_point &&
					point(to_vertex_handle(hhe_it.handle())) == rec_ec_to_point)
					break;
			}
			if (point(from_vertex_handle(hhe_it.handle())) != rec_ec_from_point ||
				point(to_vertex_handle(hhe_it.handle())) != rec_ec_to_point)
				return 0;

			// halfFaces opposite to the endpoints of the contracting edge
			HalfFaceTH hff = _hfc[from_vertex_handle(hhe_it.handle())];
			HalfFaceTH hft = _hfc[to_vertex_handle(hhe_it.handle())];
			// opposite halfFaces of hff and hft
			HalfFaceHandle ophff = hff.opposite_face_handle();
			HalfFaceHandle ophft = hft.opposite_face_handle();
			// recover opposite face relation
			if (ophff.idx() != -1)
				handle_to_entity(ophff).set_opposite_face_handle(hff.handle());
			if (ophft.idx() != -1)
				handle_to_entity(ophft).set_opposite_face_handle(hft.handle());
			// update the other two halfFaces opposite halfFace
			for (int j = 0; j < 4; j++)
			{
				if (hff.handle().idx()%4 == j || hft.handle().idx()%4 == j)
					continue;

				HalfEdgeTH he_e = handle_to_entity(_hfc[rec_ec_erase_tetra[i]*4+j].first_half_edge_handle());
				pointHandle[0] = start_point_handle(he_e.handle());
				pointHandle[1] = start_point_handle(he_e.next_half_edge_handle());
				pointHandle[2] = start_point_handle(he_e.prev_half_edge_handle());
				std::sort(pointHandle, pointHandle+3);
				sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
				if (_ophfc.find(buf) != _ophfc.end())
				{
					if (_ophfc[buf] != -1)
					{
						int idxOpHf = _ophfc[buf];
						_hfc[rec_ec_erase_tetra[i]*4+j].set_opposite_face_handle(HalfFaceHandle(idxOpHf));
						_hfc[idxOpHf].set_opposite_face_handle(HalfFaceHandle(rec_ec_erase_tetra[i]*4+j));
						_ophfc[buf] = -1;
					}
					else
					{
						std::cerr << "Error! There is more than 1 half face opposite to a particular halfface." << std::endl;
						return false;
					}
				}
				else
					_ophfc[buf] = rec_ec_erase_tetra[i]*4+j;
			}
		}
		return 1;
	}

//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------
