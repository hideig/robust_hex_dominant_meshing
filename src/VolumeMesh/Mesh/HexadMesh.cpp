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

#include <VolumeMesh/Mesh/HexadMesh.h>
#include <algorithm>
#include <set>

namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	//using namespace Hexad;

	/** Retrieve Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad to retrieve, must be valid. 
	*   \return Hexahedron which indicated by the handle
	*/
	Hexahedron& HexadMesh::handle_to_entity(HexadHandle &handle)
	{
		if ( (handle.idx() < (int)_hexac.size()) && (handle.idx() >= 0) )
		{
			return _hexac[handle.idx()];
		}
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve HalfFaceHH stored in the HexadMesh.
	*	\param  handle the index of the HalfFace to retrieve, must be valid. 
	*   \return HalfFace which indicated by the handle
	*/
	HalfFaceHH & HexadMesh::handle_to_entity(HalfFaceHandle &handle)
	{
		if ( (handle.idx() < (int)_hfc.size()) && (handle.idx() >= 0) )
		{
			return _hfc[handle.idx()];
		} 
		else
		{
			return *((HalfFaceHH*)NULL);
		}
	}

	/** Retrieve HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge to retrieve, must be valid. 
	*   \return HalfEdge which indicated by the handle
	*/
	HalfEdgeHH & HexadMesh::handle_to_entity(HalfEdgeHandle &handle)
	{
		if ( (handle.idx()< (int)_hec.size()) && (handle.idx() >= 0) )
		{
			return _hec[handle.idx()];
		} 
		else
		{
			return *((HalfEdgeHH*)NULL);
		}
	}

	/** Retrieve Vertex stored in the HexadMesh.
	*	\param  handle the index of the Vertex to retrieve, must be valid. 
	*   \return Vertex which indicated by the handle
	*/
	Vertex & HexadMesh::handle_to_entity(VertexHandle& handle)
	{
		if ( (handle.idx() < (int)_vc.size()) && (handle.idx() >= 0) )
		{
			return _vc[handle.idx()];
		} 
		else
		{
			return *((Vertex*)NULL);
		}
	}

	/** Retrieve Point stored in the HexadMesh.
	*	\param  handle the index of the Point to retrieve, must be valid. 
	*   \return Point which indicated by the handle
	*/
	Point & HexadMesh::handle_to_entity(PointHandle& handle)
	{
		if ( (handle.idx() < (int)_pc.size()) && (handle.idx() >= 0) )
		{
			return _pc[handle.idx()];
		} 
		else
		{
			return *((Point*)NULL);
		}
	}

	/** Retrieve next Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return Hexahedron which is the next of the handle.
	*/
	Hexahedron & HexadMesh::next_hexahedron(HexadHandle & handle)
	{
		if ( (handle.idx() < _npolyhedron) && (handle.idx() >= 0) )
			return handle_to_entity(next_valid_hedron_handle(handle));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve next HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFace which is the next of the handle.
	*/
	HalfFaceHH & HexadMesh::next_half_face(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle_to_entity(handle).next_half_face_handle());
	}

	/** Retrieve next HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdge which is the next of the handle.
	*/
	HalfEdgeHH & HexadMesh::next_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).next_half_edge_handle());
	}

	/** Retrieve next Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return Hexahedron which is the next of the handle.
	*/
	Hexahedron & HexadMesh::next_hexahedron(Hexahedron & hh)
	{
		if ( (hh.handle().idx() < _npolyhedron) && (hh.handle().idx() >= 0) )
			return handle_to_entity(next_valid_hedron_handle(hh.handle()));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve next HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFace which is the next of the handle.
	*/
	HalfFaceHH & HexadMesh::next_half_face(HalfFaceHH & hf)
	{
		return handle_to_entity(hf.next_half_face_handle());
	}

	/** Retrieve next HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdge which is the next of the handle.
	*/
	HalfEdgeHH & HexadMesh::next_half_edge(HalfEdgeHH & he)
	{
		return handle_to_entity(he.next_half_edge_handle());
	}

	/** Retrieve previous Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return Hexahedron which is the previous of the handle.
	*/
	Hexahedron& HexadMesh::prev_hexahedron(HexadHandle& handle)
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return handle_to_entity(prev_valid_hedron_handle(handle));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve previous HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFace which is the previous of the handle.
	*/
	HalfFaceHH& HexadMesh::prev_half_face(HalfFaceHandle & handle)
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_face_handle());
	}

	/** Retrieve previous HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdge which is the previous of the handle.
	*/
	HalfEdgeHH& HexadMesh::prev_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_edge_handle());
	}

	/** Retrieve previous Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return Hexahedron which is the previous of the handle.
	*/
	Hexahedron& HexadMesh::prev_hexahedron(Hexahedron &hh)
	{
		if ( (hh.handle().idx() <= _npolyhedron - 1)&& (hh.handle().idx() >= 0) )
			return handle_to_entity(prev_valid_hedron_handle(hh.handle()));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve previous HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFace which is the previous of the handle.
	*/
	HalfFaceHH& HexadMesh::prev_half_face(HalfFaceHH &hf)
	{
		return handle_to_entity(hf.prev_half_face_handle());
	}

	/** Retrieve previous HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdge which is the previous of the handle.
	*/
	HalfEdgeHH& HexadMesh::prev_half_edge(HalfEdgeHH& he)
	{
		return handle_to_entity(he.prev_half_edge_handle());
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad to retrieve, must be valid. 
	*   \return Hexahedron which indicated by the handle
	*/
	const Hexahedron& HexadMesh::handle_to_entity(HexadHandle &handle) const
	{
		if ( (handle.idx() < (int)_hexac.size()) && (handle.idx() >= 0) )
		{
			return _hexac[handle.idx()];
		}
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve HalfFaceHH stored in the HexadMesh.
	*	\param  handle the index of the HalfFace to retrieve, must be valid. 
	*   \return HalfFace which indicated by the handle
	*/
	const HalfFaceHH & HexadMesh::handle_to_entity(HalfFaceHandle &handle) const
	{
		if ( (handle.idx() < (int)_hfc.size()) && (handle.idx() >= 0) )
		{
			return _hfc[handle.idx()];
		} 
		else
		{
			return *((HalfFaceHH*)NULL);
		}
	}

	/** Retrieve HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge to retrieve, must be valid. 
	*   \return HalfEdge which indicated by the handle
	*/
	const HalfEdgeHH & HexadMesh::handle_to_entity(HalfEdgeHandle &handle) const
	{
		if ( (handle.idx()< (int)_hec.size()) && (handle.idx() >= 0) )
		{
			return _hec[handle.idx()];
		} 
		else
		{
			return *((HalfEdgeHH*)NULL);
		}
	}

	/** Retrieve Vertex stored in the HexadMesh.
	*	\param  handle the index of the Vertex to retrieve, must be valid. 
	*   \return Vertex which indicated by the handle
	*/
	const Vertex & HexadMesh::handle_to_entity(VertexHandle& handle) const
	{
		if ( (handle.idx() < (int)_vc.size()) && (handle.idx() >= 0) )
		{
			return _vc[handle.idx()];
		} 
		else
		{
			return *((Vertex*)NULL);
		}
	}

	/** Retrieve Point stored in the HexadMesh.
	*	\param  handle the index of the Point to retrieve, must be valid. 
	*   \return Point which indicated by the handle
	*/
	const Point & HexadMesh::handle_to_entity(PointHandle& handle) const
	{
		if ( (handle.idx() < (int)_pc.size()) && (handle.idx() >= 0) )
		{
			return _pc[handle.idx()];
		} 
		else
		{
			return *((Point*)NULL);
		}
	}

	/** Retrieve next Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return Hexahedron which is the next of the handle.
	*/
	const Hexahedron & HexadMesh::next_hexahedron(HexadHandle & handle) const
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return handle_to_entity(next_valid_hedron_handle(handle));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve next HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFace which is the next of the handle.
	*/
	const HalfFaceHH & HexadMesh::next_half_face(HalfFaceHandle& handle) const
	{
		return handle_to_entity(handle_to_entity(handle).next_half_face_handle());
	}

	/** Retrieve next HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdge which is the next of the handle.
	*/
	const HalfEdgeHH & HexadMesh::next_half_edge(HalfEdgeHandle &handle) const
	{
		return handle_to_entity(handle_to_entity(handle).next_half_edge_handle());
	}

	/** Retrieve next Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return Hexahedron which is the next of the handle.
	*/
	const Hexahedron & HexadMesh::next_hexahedron(Hexahedron & hh) const
	{
		if ( (hh.handle().idx() <= _npolyhedron - 1) && (hh.handle().idx() >= 0) )
			return handle_to_entity(next_valid_hedron_handle(hh.handle()));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve next HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFace which is the next of the handle.
	*/
	const HalfFaceHH & HexadMesh::next_half_face(HalfFaceHH & hf) const
	{
		return handle_to_entity(hf.next_half_face_handle());
	}

	/** Retrieve next HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdge which is the next of the handle.
	*/
	const HalfEdgeHH & HexadMesh::next_half_edge(HalfEdgeHH & he) const
	{
		return handle_to_entity(he.next_half_edge_handle());
	}

	/** Retrieve previous Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return Hexahedron which is the previous of the handle.
	*/
	const Hexahedron& HexadMesh::prev_hexahedron(HexadHandle& handle) const
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return handle_to_entity(prev_valid_hedron_handle(handle));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve previous HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFace which is the previous of the handle.
	*/
	const HalfFaceHH& HexadMesh::prev_half_face(HalfFaceHandle & handle) const
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_face_handle());
	}

	/** Retrieve previous HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdge which is the previous of the handle.
	*/
	const HalfEdgeHH& HexadMesh::prev_half_edge(HalfEdgeHandle &handle) const
	{
		return handle_to_entity(handle_to_entity(handle).prev_half_edge_handle());
	}

	/** Retrieve previous Hexad stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return Hexahedron which is the previous of the handle.
	*/
	const Hexahedron& HexadMesh::prev_hexahedron(Hexahedron &hh) const
	{
		if ( (hh.handle().idx() <= _npolyhedron - 1)&& (hh.handle().idx() >= 0) )
			return handle_to_entity(prev_valid_hedron_handle(hh.handle()));
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve previous HalfFace stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFace which is the previous of the handle.
	*/
	const HalfFaceHH& HexadMesh::prev_half_face(HalfFaceHH &hf) const
	{
		return handle_to_entity(hf.prev_half_face_handle());
	}

	/** Retrieve previous HalfEdge stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdge which is the previous of the handle.
	*/
	const HalfEdgeHH& HexadMesh::prev_half_edge(HalfEdgeHH& he) const
	{
		return handle_to_entity(he.prev_half_edge_handle());
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve next HexadHandle stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return HexadHandle which is the next hexahedron handle.
	*/
	HexadHandle  HexadMesh::next_hexahedron_handle(HexadHandle & handle)
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return next_valid_hedron_handle(handle);
		else
			return (HexadHandle(-1));
	}

	/** Retrieve next HalfFaceHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFaceHandle which is the next half face handle.
	*/
	HalfFaceHandle  HexadMesh::next_half_face_handle(HalfFaceHandle& handle)
	{
		return handle_to_entity(handle).next_half_face_handle();
	}

	/** Retrieve next HalfEdgeHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdgeHandle which is the next half edge handle.
	*/
	HalfEdgeHandle  HexadMesh::next_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).next_half_edge_handle();
	}

	/** Retrieve next HexadHandle stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its next, must be valid. 
	*   \return HexadHandle which is the next hexahedron handle.
	*/
	HexadHandle  HexadMesh::next_hexahedron_handle(Hexahedron & hh)
	{
		if ( (hh.handle().idx() <= _npolyhedron - 1) && (hh.handle().idx() >= 0) )
			return next_valid_hedron_handle(hh.handle());
		else
			return HexadHandle(-1);
	}

	/** Retrieve next HalfFaceHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its next, must be valid. 
	*   \return HalfFaceHandle which is the next half face handle.
	*/
	HalfFaceHandle  HexadMesh::next_half_face_handle(HalfFaceHH & hf)
	{
		return hf.next_half_face_handle();
	}

	/** Retrieve next HalfEdgeHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its next, must be valid. 
	*   \return HalfEdgeHandle which is the next half edge handle.
	*/
	HalfEdgeHandle  HexadMesh::next_half_edge_handle(HalfEdgeHH & he)
	{
		return he.next_half_edge_handle();
	}

	/** Retrieve previous HexadHandle stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return HexadHandle which is the previous hexahedron handle.
	*/
	HexadHandle HexadMesh::prev_hexahedron_handle(HexadHandle& handle)
	{
		if ( (handle.idx() <= _npolyhedron - 1) && (handle.idx() >= 0) )
			return prev_valid_hedron_handle(handle);
		else
			return HexadHandle(-1);
	}

	/** Retrieve previous HalfFaceHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFaceHandle which is the previous half face handle.
	*/
	HalfFaceHandle HexadMesh::prev_half_face_handle(HalfFaceHandle & handle)
	{
		return handle_to_entity(handle).prev_half_face_handle();
	}

	/** Retrieve previous HalfEdgeHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdge which is the previous half edge handle.
	*/
	HalfEdgeHandle HexadMesh::prev_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).prev_half_edge_handle();
	}

	/** Retrieve previous Hexad handle stored in the HexadMesh.
	*	\param  handle the index of the Hexad which we want to get its previous, must be valid. 
	*   \return HexadHandle which is the previous hexahedron handle.
	*/
	HexadHandle HexadMesh::prev_hexahedron_handle(Hexahedron &hh)
	{
		if ( (hh.handle().idx() <= _npolyhedron - 1)&& (hh.handle().idx() >= 0) )
			return prev_valid_hedron_handle(hh.handle());
		else
			return HexadHandle(-1);
	}

	/** Retrieve previous HalfFaceHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfFace which we want to get its previous, must be valid. 
	*   \return HalfFaceHandle which is the previous half face handle.
	*/
	HalfFaceHandle HexadMesh::prev_half_face_handle(HalfFaceHH &hf)
	{
		return hf.prev_half_face_handle();
	}

	/** Retrieve previous HalfEdgeHandle stored in the HexadMesh.
	*	\param  handle the index of the HalfEdge which we want to get its previous, must be valid. 
	*   \return HalfEdgeHandle which is the previous half edge handle.
	*/
	HalfEdgeHandle HexadMesh::prev_half_edge_handle(HalfEdgeHH& he)
	{
		return he.prev_half_edge_handle();
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve the first half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::first_hedron_half_edge_handle(HedronIter & hIter_)
	{
		assert(hIter_.handle().is_valid());
		return HalfEdgeHandle(hIter_.handle().idx() * 12);
	}
	/** Retrieve the first half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::first_hedron_half_edge_handle(HedronHandle & hh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12);
	}

	extern const int HMate[6][8];
	
	/** Retrieve the first half face handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::first_vertex_half_face_handle(VertexHandle & vh_)
	{
		assert(vh_.is_valid());
		int idx;
		for (int i = 0; i < 6; i++)
		{
			if (HMate[i][vh_%8] != -1)
			{
				idx = HMate[i][vh_%8];
				break;
			}
		}
		return HalfFaceHandle((vh_>>3)*6+idx);
	}
		
	/** Retrieve the first half face handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::first_vertex_half_face_handle(VertexIter   & vIter_)
	{
		VertexHandle vh_ = vIter_.handle();
		assert(vh_.is_valid());
		int idx;
		for (int i = 0; i < 6; i++)
		{
			if (HMate[i][vh_%8] != -1)
			{
				idx = HMate[i][vh_%8];
				break;
			}
		}
		return HalfFaceHandle((vh_>>3)*6+idx);
	}

	/** Retrieve the first half edge handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::first_vertex_half_edge_handle(VertexHandle & vh_)
	{
		assert(vh_.is_valid());
		return HalfEdgeHandle((vh_>>3)*24+vhecHH_[vh_%8]);
	}

	/** Retrieve the first half edge handle related to a vertex
	*   \param in vh_ the vertex handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::first_vertex_half_edge_handle(VertexIter   & vIter_)
	{
		assert(vIter_.handle().is_valid());
		return HalfEdgeHandle((vIter_.handle()>>3)*24+vhecHH_[vIter_.handle()%8]);
	}

	/** Retrieve the next half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::next_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_)
	{
		assert(hIter_.handle().is_valid() && heh_.is_valid());
		return HalfEdgeHandle(hIter_.handle() * 12 + (heh_.idx() + 1) % 12);
	}

	/** Retrieve the next half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::next_hedron_half_edge_handle(HedronHandle & hh_, HalfEdgeHandle & heh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12 + (heh_.idx() + 1) % 12);
	}

	/** Retrieve the previous half edge handle from a hedron iterator
	*   \param in hIter_ the iterator
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::prev_hedron_half_edge_handle(HedronIter & hIter_, HalfEdgeHandle & heh_)
	{
		assert(hIter_.handle().is_valid() && heh_.is_valid());
		return HalfEdgeHandle(hIter_.handle() * 12 + (heh_.idx() - 1) % 12);
	}

	/** Retrieve the previous half edge handle from a hedron iterator
	*   \param in hh_ the hedron handle
	*   \param in heh_ the current half edge handle
	*   \return the corresponding half edge handle
	*/
	HalfEdgeHandle HexadMesh::prev_hedron_half_edge_handle(HedronHandle & hh_, HalfEdgeHandle & heh_)
	{
		assert(hh_.is_valid());
		return HalfEdgeHandle(hh_.idx() * 12 + (heh_.idx() - 1) % 12);
	}

	/** Retrieve the next half face handle relatec to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::next_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_)
	{
		assert(vh_.is_valid());
		return HalfFaceHandle((vh_>>3)*6+HMate[hfh_%6][vh_%8]);
	}

	/** Retrieve the next half face handle relatec to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::next_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_)
	{
		assert(vIter_.handle().is_valid());
		return HalfFaceHandle((vIter_.handle()>>3)*6+HMate[hfh_%6][vIter_.handle()%8]);
	}

	/** Retrieve the previous half face handle relatec to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::prev_vertex_half_face_handle(VertexHandle & vh_,  HalfFaceHandle & hfh_)
	{
		assert(vh_.is_valid());
		return HalfFaceHandle((vh_>>3)*6+HMate[HMate[hfh_%6][vh_%8]][vh_%8]);
	}
	
	/** Retrieve the previous half face handle relatec to a vertex
	*   \param in vh_ the vertex handle
	*   \param in hfh_ the current half face handle
	*   \return the corresponding half face handle
	*/
	HalfFaceHandle HexadMesh::prev_vertex_half_face_handle(VertexIter & vIter_, HalfFaceHandle & hfh_)
	{
		assert(vIter_.handle().is_valid());
		return HalfFaceHandle((vIter_.handle()>>3)*6+HMate[HMate[hfh_%6][vIter_.handle()%8]][vIter_.handle()%8]);
	}
	//-----------------------------------------------------------------------------------------------------------------


	/** Access from HalfEdge to its Point.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return PointHandle which gives the start point of HalfEdge.
	*/
	PointHandle HexadMesh::start_point_handle(HalfEdgeHandle &handle)
	{
		if ( (handle.idx() < (int)_hec.size()) && (handle.idx() >= 0) )
		{			
			return _vc[_hec[handle.idx()].start_vertex_handle().idx()].point_handle();
		} 
		else
		{
			return PointHandle(-1);
		}
	}

	/** Access from HalfEdge to its Point.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return Point which is the start point of HalfEdge.
	*/
	Point& HexadMesh::start_point(HalfEdgeHandle &handle)
	{
		return handle_to_entity(start_point_handle(handle));
	}

	/** Get Radial HalfEdge of a particular HalfEdge.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return Radial HalfEdge.
	*/
	HalfEdgeHH & HexadMesh::radial_half_edge(HalfEdgeHandle &handle)
	{
		if ( (handle.idx() < (int)_hec.size()) && (handle.idx() >= 0) )
		{
			HalfEdgeHH he = handle_to_entity(handle);
			HalfFaceHH hf  = handle_to_entity(he.half_face_handle());
			if (hf.has_opposite_face())
			{
				PointHandle nextp = start_point_handle(he.next_half_edge_handle());
				HalfFaceHH opphf = opposite_half_face(hf.handle());
				HalfEdgeHH opphe = handle_to_entity(opphf.first_half_edge_handle());
				while( nextp.idx() != start_point_handle(opphe.handle()).idx() )
				{
					opphe = handle_to_entity(opphe.next_half_edge_handle());
				}
				return handle_to_entity(opphe.handle());
			}
			else
				return *((HalfEdgeHH*)NULL);
		} 
		else
		{
			return *((HalfEdgeHH*)NULL);
		}
	}

	/** Get Radial HalfEdgeHandle of a particular HalfEdge.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return Radial HalfEdgeHandle.
	*/
	HalfEdgeHandle  HexadMesh::radial_half_edge_handle(HalfEdgeHandle &handle)
	{
		if ( (handle.idx() < (int)_hec.size()) && (handle.idx() >= 0) )
		{
			HalfEdgeHH he = handle_to_entity(handle);
			HalfFaceHH hf  = handle_to_entity(he.half_face_handle());
			if (hf.has_opposite_face())
			{
				PointHandle nextp = start_point_handle(he.next_half_edge_handle());
				HalfFaceHH opphf = opposite_half_face(hf.handle());
				HalfEdgeHH opphe = handle_to_entity(opphf.first_half_edge_handle());
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

	/** Check whether the HalfEdge has its Radial HalfEdge.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::has_radial_half_edge(HalfEdgeHandle &handle)
	{
		HalfEdgeHH he = handle_to_entity(handle);
		HalfFaceHH hf  = handle_to_entity(he.half_face_handle());
		return hf.has_opposite_face();
	}

	/** Get mate HalfEdge of a particular HalfEdge.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return Mate HalfEdge.
	*/
	HalfEdgeHH& HexadMesh::mate_half_edge(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle_to_entity(handle).mate_half_edge_handle());
	}

	/** Get mate HalfEdgeHandle of a particular HalfEdge.
	*	\param  handle the index of the HalfEdge, must be valid. 
	*   \return Mate HalfEdgeHandle.
	*/
	HalfEdgeHandle HexadMesh::mate_half_edge_handle(HalfEdgeHandle &handle)
	{
		return handle_to_entity(handle).mate_half_edge_handle();
	}

	/** Get Opposite HalfFace of a particular HalfFace.
	*	\param  handle the index of the HalfFace, must be valid. 
	*   \return Opposite HalfFace.
	*/
	HalfFaceHH& HexadMesh::opposite_half_face(HalfFaceHandle& handle)
	{
		if ( (handle.idx() < (int)_hec.size()) && (handle.idx() >= 0) )
		{
			HalfFaceHH hf = handle_to_entity(handle);
			if (hf.has_opposite_face())
				return handle_to_entity(hf.opposite_face_handle());
			else
				return *((HalfFaceHH*)NULL);
		}
		else
		{
			return *((HalfFaceHH*)NULL);
		}
	}

	/** Get Opposite HalfFace of a particular HalfFace.
	*	\param  handle the index of the HalfFace, must be valid. 
	*   \return Opposite HalfFace.
	*/
	HalfFaceHandle HexadMesh::opposite_half_face_handle(HalfFaceHandle& handle)
	{
		if ( (handle.idx() < (int)_hec.size()) && (handle.idx() >= 0) )
		{
			HalfFaceHH hf = handle_to_entity(handle);
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

	/** Check whether the HalfFace has its Opposite HalfFace.
	*	\param  handle the index of the HalfFace, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::has_opposite_half_face(HalfFaceHandle & handle)
	{
		return handle_to_entity(handle).has_opposite_face();
	}

	/** Get Hexahedron of a particular Vertex.
	*	\param  handle the index of the Vertex, must be valid. 
	*   \return Hexahedron.
	*/
	Hexahedron& HexadMesh::hexahedron(VertexHandle& handle)
	{
		return handle_to_entity(HexadHandle(handle.idx() / 8));
	}

	/** Get half edge handle with a particular face and from vertex.
	*	\param  fh the index of the half face, must be valid. 
	*	\param  vh the index of the vertex, must be valid. 
	*   \return half edge handle.
	*/
	HalfEdgeHandle & HexadMesh::half_edge_handle(HalfFaceHandle & fh, VertexHandle & vh)
	{
		HalfEdgeHH he = handle_to_entity(handle_to_entity(fh).first_half_edge_handle());
		for ( int i = 0; i < 4; i ++ )
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
	HalfFaceHandle & HexadMesh::half_face_handle(HalfEdgeHandle & handle)
	{
		return handle_to_entity(handle).half_face_handle();
	}

	/** Check whether the HalfFace has particular Vertex.
	*	\param  hfHandle the index of the HalfFace, must be valid. 
	*	\param  vHandle the index of the Vertex, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::has_vertex(HalfFaceHandle hfHandle, VertexHandle vHandle)
	{
		HalfFaceHH hf = handle_to_entity(hfHandle);
		HalfEdgeHH he = handle_to_entity(hf.first_half_edge_handle());
		for ( int i = 0; i < 4; i ++ )
		{
			if ( he.start_vertex_handle() == vHandle )
				return true;
			he = next_half_edge(he.handle());
		}
		return false;
	}

	/** Check whether the HalfFace has particular Point.
	*	\param  hfHandle the index of the HalfFace, must be valid. 
	*	\param  pHandle the index of the Point, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::has_point(HalfFaceHandle hfHandle, PointHandle pHandle)
	{
		HalfFaceHH hf = handle_to_entity(hfHandle);
		HalfEdgeHH he = handle_to_entity(hf.first_half_edge_handle());
		for ( int i = 0; i < 4; i++ )
		{
			if ( start_point_handle(he.handle()) == pHandle )
				return true;
			he = next_half_edge(he.handle());
		}
		return false;
	}


	/** Check whether the hedron handle has boundary half face.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::has_boundary_face(PolyhedronHandle & hh_)
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
	bool HexadMesh::has_boundary_face(HedronIter & hIter_)
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




	/** Release all the topology stored in the HexadMesh class.
	*/
	void HexadMesh::release_mesh()
	{
		_hec.clear();
		_hfc.clear();
		_hexac.clear();
		_pc.clear();
		_vc.clear();
		_ophfc.clear();
		_npolyhedron = 0;
		_npoint = 0;
	}
		
	/** Add a point to the TetrasMesh.
	*   \return new point's PointHandle
	*/
	PointHandle HexadMesh::add_point(Point &p)
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
	int HexadMesh::erase_point(PointHandle &handle)
	{
		//invalid handle
		if ( handle.idx() < 0 || handle.idx() >= _npoint || !is_valid(handle) )
			return -1;

		// add pointHandle to free_point_list.
		add_to_free_point_list(handle);

		//--------------------------new added---------------------------------------//
		//_pvc[ph.idx()] = VertexHandle(-1);
		//--------------------------end added---------------------------------------//

		return 1;
	}

	/** Add a Hexad to the HexadsMesh.
	*	\param  pntVec point of the hexad to add. 
	*   \return -1 invalid parameter.
	*			 0 already has the hexad.
	*			 1 success.
	*/
	HexadHandle HexadMesh::add_hedron(std::vector<Point> & pntVec)
	{
		if ( pntVec.size()!= 8 )
			return HexadHandle(-1);

		int i;	//this is insert hexahedron index
		int hfIdx;
		int heIdx;
		int verIdx;
		int pidx;
		int _npoint_old = _npoint;
		HexadHandle hh;
		std::vector<Point>::iterator iter;
		// insert into empty location, just pop free_list and adjust pointHandle
		if ( _free_hedron_list.empty()==false )
		{
			hh = pop_erase_from_free_tetra_list();
			i = hh.idx();
			hfIdx =i * 6;
			heIdx = i * 24;
			verIdx = i * 8;

			//reset halfFace again to avoid invalid case
			for ( int j = 0; j < 6; j ++ )
				_hfc[hfIdx+j] = HalfFaceHH( i * 6 + j);
			// Halfedge of halfFace 0
			_hec[heIdx]    = HalfEdgeHH(heIdx,   verIdx+2);
			_hec[heIdx+1]  = HalfEdgeHH(heIdx+1, verIdx+3);
			_hec[heIdx+2]  = HalfEdgeHH(heIdx+2, verIdx+7);
			_hec[heIdx+3]  = HalfEdgeHH(heIdx+3, verIdx+6);
			// Halfedge of halfFace 1
			_hec[heIdx+4]  = HalfEdgeHH(heIdx+4, verIdx+4);
			_hec[heIdx+5]  = HalfEdgeHH(heIdx+5, verIdx+5);
			_hec[heIdx+6]  = HalfEdgeHH(heIdx+6, verIdx+6);
			_hec[heIdx+7]  = HalfEdgeHH(heIdx+7, verIdx+7);
			// Halfedge of halfFace 2
			_hec[heIdx+8]  = HalfEdgeHH(heIdx+8, verIdx+1);
			_hec[heIdx+9]  = HalfEdgeHH(heIdx+9, verIdx+2);
			_hec[heIdx+10] = HalfEdgeHH(heIdx+10, verIdx+6);
			_hec[heIdx+11] = HalfEdgeHH(heIdx+11, verIdx+5);
			// Halfedge of halfFace 3
			_hec[heIdx+12] = HalfEdgeHH(heIdx+12, verIdx+0);
			_hec[heIdx+13] = HalfEdgeHH(heIdx+13, verIdx+3);
			_hec[heIdx+14] = HalfEdgeHH(heIdx+14, verIdx+2);
			_hec[heIdx+15] = HalfEdgeHH(heIdx+15, verIdx+1);
			// Halfedge of halfFace 4
			_hec[heIdx+16] = HalfEdgeHH(heIdx+16, verIdx+0);
			_hec[heIdx+17] = HalfEdgeHH(heIdx+17, verIdx+4);
			_hec[heIdx+18] = HalfEdgeHH(heIdx+18, verIdx+7);
			_hec[heIdx+19] = HalfEdgeHH(heIdx+19, verIdx+3);
			// Halfedge of halfFace 5
			_hec[heIdx+20] = HalfEdgeHH(heIdx+20, verIdx+0);
			_hec[heIdx+21] = HalfEdgeHH(heIdx+21, verIdx+1);
			_hec[heIdx+22] = HalfEdgeHH(heIdx+22, verIdx+5);
			_hec[heIdx+23] = HalfEdgeHH(heIdx+23, verIdx+4);
		}
		// insert at tail
		else
		{
			i = _hexac.size();
			hfIdx =i * 6;
			heIdx = i * 24;
			verIdx = i * 8;

			_hexac.push_back(Hexahedron(i));
			for ( int j=0; j<6; j++ )
				_hfc.push_back(HalfFaceHH(hfIdx + j));
			for ( int j=0; j<8; j++ )
				_vc.push_back(Vertex(-1));
			// Halfedge of halfface 0
			_hec.push_back(HalfEdgeHH(heIdx,   verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+1, verIdx+3));
			_hec.push_back(HalfEdgeHH(heIdx+2, verIdx+7));
			_hec.push_back(HalfEdgeHH(heIdx+3, verIdx+6));
			// Halfedge of halfface 1
			_hec.push_back(HalfEdgeHH(heIdx+4, verIdx+4));
			_hec.push_back(HalfEdgeHH(heIdx+5, verIdx+5));
			_hec.push_back(HalfEdgeHH(heIdx+6, verIdx+6));
			_hec.push_back(HalfEdgeHH(heIdx+7, verIdx+7));
			// Halfedge of halfface 2
			_hec.push_back(HalfEdgeHH(heIdx+8, verIdx+1));
			_hec.push_back(HalfEdgeHH(heIdx+9, verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+10, verIdx+6));
			_hec.push_back(HalfEdgeHH(heIdx+11, verIdx+5));
			// Halfedge of halfface 3
			_hec.push_back(HalfEdgeHH(heIdx+12, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+13, verIdx+3));
			_hec.push_back(HalfEdgeHH(heIdx+14, verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+15, verIdx+1));
			// Halfedge of halfface 4
			_hec.push_back(HalfEdgeHH(heIdx+16, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+17, verIdx+4));
			_hec.push_back(HalfEdgeHH(heIdx+18, verIdx+7));
			_hec.push_back(HalfEdgeHH(heIdx+19, verIdx+3));
			// Halfedge of halfface 5
			_hec.push_back(HalfEdgeHH(heIdx+20, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+21, verIdx+1));
			_hec.push_back(HalfEdgeHH(heIdx+22, verIdx+5));
			_hec.push_back(HalfEdgeHH(heIdx+23, verIdx+4));
		}

		//find & add the vertex & point first
		PointFindIF<Point> pif(pntVec[0]);
		for (unsigned int j = 0; j < pntVec.size(); j ++)
		{
			pif = pntVec[j];
			// add new point
			if( (pidx=find_point(pif).idx())==-1 )
			{
				//add at tail
				if ( _free_point_list.empty()==true )
				{
					_pvc.push_back(VertexHandle(i*8+j));
					_pc.push_back(pntVec[j]);
					_vc[i*8+j] = Vertex(i*8+j, _npoint);
					//_pvm[PointHandle(_npoint)].push_back(VertexHandle(i*8+j));
					this->_npoint++;
				}
				//add at free list location
				else
				{
					PointHandle ph = pop_erase_from_free_point_list();
					Point& p = handle_to_entity(ph);
					p = pntVec[j];
					_vc[i*8+j] = Vertex(i*8+j, ph.idx());
					_pvc[ph] = VertexHandle(i*4+3);
					//_pvm[ph].push_back(VertexHandle(i*8+j));
				}
			}
			// old point
			else
			{
				_vc[i*8+j] = Vertex(i*8, pidx);
				_vc[i*8+j].set_next_vertex_handle(_pvc[pidx]);
				_pvc[pidx] = _vc[i*8+j].handle();
				//_pvm[PointHandle(pidx)].push_back(VertexHandle(i*8+j));
			}
		}

		int pointHandle[4];
		char buf[6][128];
		int ophfRollBack[6];
		for ( int j=0; j<6; j++ )
		{
			ophfRollBack[j] = -1;
		}

		//dealing with half face shared with old tetra
		for ( int j = 0; j < 6; j++ )
		{
			pointHandle[0] = start_point_handle(HalfEdgeHandle(heIdx + j*4)).idx();
			pointHandle[1] = start_point_handle(HalfEdgeHandle(heIdx + j*4 + 1)).idx();
			pointHandle[2] = start_point_handle(HalfEdgeHandle(heIdx + j*4 + 2)).idx();
			pointHandle[3] = start_point_handle(HalfEdgeHandle(heIdx + j*4 + 3)).idx();
			std::sort(pointHandle, pointHandle + 4);
			sprintf_s(buf[j],"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
			//find its opposite half face
			if ( _ophfc.find(buf[j])!=_ophfc.end() )
			{
				if ( _ophfc[buf[j]]!=-1 )
				{
					int idxOpHf = _ophfc[buf[j]];
					_hfc[idxOpHf].set_opposite_face_handle(HalfFaceHandle(hfIdx+j));
					_hfc[hfIdx+j].set_opposite_face_handle(HalfFaceHandle(idxOpHf));
					_ophfc[buf[j]] = -1;
					ophfRollBack[j] = idxOpHf;
				} 
				else
				{
					std::cerr<<"Error: There is more than 1 half face opposite to another."<< std::endl;
					if (i==_hexac.size()-1)
					{
						_hexac.pop_back();
						for ( int k=0; k<6; k++ )
							_hfc.pop_back();
						for ( int k=0; k<24; k++ )
							_hec.pop_back();					
						for ( int k=0; k<8; k++ )
						{
							_pvc[_vc[i*8+k].point_handle().idx()] = _vc[i*8+k].next_vertex_handle();
							_vc.pop_back();
						}
					}
					else
					{					
						for ( int k=0; k<8; k++ )
							_pvc[_vc[i*8+k].point_handle().idx()] = _vc[i*8+k].next_vertex_handle();
						add_to_free_hedron_list(_hexac[i].handle());
					}

					// _pc roll back
					for (int k=0; k<_npoint-_npoint_old; k++)
					{
						_pc.pop_back();
						_pvc.pop_back();
					}

					//here we introduce _ophfc roll back
					//such function is not support in tetraMesh now
					for ( int k = 0; k < j; k++ )
					{
						//need roll back for _ophfc
						if ( ophfRollBack[k]!=-1 )
						{
							_ophfc[buf[k]] = ophfRollBack[k];
							_hfc[ophfRollBack[k]].set_opposite_face_handle(HalfFaceHandle(-1));
						}
						else
							_ophfc.erase(buf[k]);
					}
					_npoint = _pc.size();
					return HexadHandle(-1);
				}
			}
			else
			{
				_ophfc[buf[j]] = hfIdx + j;
			}
		}
		this->_npolyhedron++;
		return HexadHandle(i);
	}

	/** Erase a Hexad of the HexadMesh.
	*	This function exchange the last hexad and the hexad to delete,
	*	then delete the last hexad and update the opposite relation.
	*	\param  handle index of the hexad to erase. 
	*   \return -1 invalid parameter.
	*			 1 success.
	*/
	HexadHandle HexadMesh::erase_hedron(VolumeMesh::HexadHandle &handle)
	{
		//using namespace std;
		//invalid handle
		if ( (handle.idx() < 0) || (handle.idx() >= _npolyhedron) )
		{
			std::cerr << "Error: HexadMesh::erase_hedron/HexadHandle is out of size!" <<std::endl;
			return HexadHandle(-1);
		}

		int hfIdx =handle * 6;
		int heIdx = handle * 24;
		int verIdx = handle * 8;

		//reset opposite halfface relation & delete point if necessary
		Hexahedron hexad = handle_to_entity(handle);
		std::set<PointHandle> ps;
		for ( int i = 0; i < 6; i ++ )
		{
			char buf[128];
			PointHandle pointHandle[4];
			HalfFaceHH hf = handle_to_entity(HalfFaceHandle(hexad.first_half_face_handle().idx() + i));
			if ( hf.has_opposite_face() )
			{
				HalfFaceHH & opphf = handle_to_entity(hf.opposite_face_handle());
				HalfEdgeHandle oppheHandle = opphf.first_half_edge_handle();
				//reset opposite halfface relation
				opphf.set_opposite_face_handle(HalfFaceHandle(-1));
				//reset opposite map
				pointHandle[0] = start_point_handle(HalfEdgeHandle(oppheHandle));
				pointHandle[1] = start_point_handle(HalfEdgeHandle(oppheHandle+1));
				pointHandle[2] = start_point_handle(HalfEdgeHandle(oppheHandle+2));
				pointHandle[3] = start_point_handle(HalfEdgeHandle(oppheHandle+3));
				std::sort(pointHandle, pointHandle+4);
				sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
				_ophfc[buf] = opphf.handle().idx();					
				ps.insert(pointHandle[0]);
				ps.insert(pointHandle[1]);
				ps.insert(pointHandle[2]);
				ps.insert(pointHandle[3]);
			}
			else
			{
				HalfEdgeHandle heHandle = hf.first_half_edge_handle();
				//erase opposite map relation
				pointHandle[0] = start_point_handle(HalfEdgeHandle(heHandle));
				pointHandle[1] = start_point_handle(HalfEdgeHandle(heHandle + 1));
				pointHandle[2] = start_point_handle(HalfEdgeHandle(heHandle + 2));
				pointHandle[3] = start_point_handle(HalfEdgeHandle(heHandle + 3));
				std::sort(pointHandle,pointHandle+4);
				sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
				_ophfc.erase(buf);
			}
		}
		//delete point
		if ( ps.size() < 8 )
		{
			for ( int i = 0; i < 8; i ++ )
			{
				Vertex v = handle_to_entity(VertexHandle(hexad.first_vertex_handle().idx() + i));
				// This point might be delete. Cuz It doesn't belong to opposite tetra.
				if ( ps.find(v.point_handle()) == ps.end() )
				{
					VertexFindIF vif(v);
					if ( find_if(_vc.begin(), _vc.begin() + hexad.first_vertex_handle().idx(), vif)!=(_vc.begin() + hexad.first_vertex_handle().idx()) )
						continue;
					if ( find_if(_vc.begin()+hexad.first_vertex_handle()+8,_vc.end(),vif)!=_vc.end() )
						continue;
					// This point has no reference, should be delete. 2 rule.
					// #1 We exchange it with the last point in vector, then pop vector.
					// #2 The erase point is the last point in vector, just pop it.
					// #3 Revise the opposite mapping relation
					if ( v.point_handle().idx() != _npoint-1 )
					{
						for (unsigned int j=0; j < _hfc.size(); j ++)
						{
							if ( has_point(HalfFaceHandle(j), PointHandle(_npoint-1)) )
							{
								char buf[128];
								int tmpHf;

								//Revise opposite mapping relation
								PointHandle pointHandle[4];
								pointHandle[0] = start_point_handle(HalfEdgeHandle(j*4));
								pointHandle[1] = start_point_handle(HalfEdgeHandle(j*4+1));
								pointHandle[2] = start_point_handle(HalfEdgeHandle(j*4+2));
								pointHandle[3] = start_point_handle(HalfEdgeHandle(j*4+3));
								std::sort(pointHandle, pointHandle+4);
								sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
								if ( _ophfc.find(buf)!=_ophfc.end() )
								{
									tmpHf = _ophfc[buf];
									_ophfc.erase(buf);

									pointHandle[3] = v.point_handle();
									std::sort(pointHandle, pointHandle+4);
									sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
									_ophfc[buf] = tmpHf;
								}
							}
						}
						add_to_free_point_list(v.point_handle());
					}

					//--------------------------new added---------------------------------------//
					_pvc[v.point_handle().idx()] = VertexHandle(-1);
					//--------------------------end added---------------------------------------//
				}
			}
		}
		//Hexad handle is not the last hexad,need to exchange hexad.
		//if ( handle.idx() != _npolyhedron - 1 )
		//{
		//	Hexahedron lastHexad = handle_to_entity(HexadHandle(_npolyhedron - 1));
		//	HalfFaceHH lastHf = handle_to_entity(lastHexad.first_half_face_handle());
		//	HalfFaceHH hf = handle_to_entity(hexad.first_half_face_handle());
		//	
		//	//reset last halfface's opposite relation to new position
		//	for ( int i = 0; i < 6; i ++ )
		//	{
		//		if ( lastHf.has_opposite_face() )
		//		{
		//			HalfFaceHH& opphf = handle_to_entity(lastHf.opposite_face_handle());
		//			opphf.set_opposite_face_handle(HalfFaceHandle(handle.idx() * 6 + i));
		//			_hfc[hf.handle()].set_opposite_face_handle(opphf.handle());
		//		}
		//		else
		//		{
		//			_hfc[hf.handle()].set_opposite_face_handle(HalfFaceHandle(-1));

		//			// reset _ophfc' Map relation for the last Hexad
		//			char buf[128];
		//			PointHandle pointHandle[4];
		//			HalfEdgeHandle lastHe = lastHf.first_half_edge_handle();
		//			pointHandle[0] = start_point_handle(HalfEdgeHandle(lastHe));
		//			pointHandle[1] = start_point_handle(HalfEdgeHandle(lastHe+1));
		//			pointHandle[2] = start_point_handle(HalfEdgeHandle(lastHe+2));
		//			pointHandle[3] = start_point_handle(HalfEdgeHandle(lastHe+3));
		//			std::sort(pointHandle, pointHandle+4);
		//			sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
		//			_ophfc[buf] = handle*6+i;
		//		}
		//		hf = next_half_face(hf);
		//		lastHf = next_half_face(lastHf);
		//	}

		//	//reset point handle 
		//	for ( int i=0; i<8; i++ )
		//	{
		//		Vertex lastV = handle_to_entity(VertexHandle(lastHexad.first_vertex_handle().idx()+i));
		//		Vertex &v = handle_to_entity(VertexHandle(hexad.first_vertex_handle().idx()+i));
		//		v.set_point_handle(lastV.point_handle());
		//	}
		//}
		//delete all the topology relation
		//_hexac.pop_back();
		//for ( int i=0; i<6; i++ )
		//	_hfc.pop_back();
		//for ( int i=0;i<24;i++ )
		//	_hec.pop_back();
		//for ( int i=0; i<8; i++ )
		//	_vc.pop_back();
		//_npolyhedron = _hexac.size();
		add_to_free_hedron_list(handle);
		return next_valid_hedron_handle(handle);
	}

	/** Retrieve EdgeStar in HexadMesh.
	*	\param  eh handle of the HalfEdge. 
	*	\param  hexadVec store the result of the EdgeStar. 
	*   \return number of the HalfEdge in EdgeStar.
	*/
	int HexadMesh::edge_star(HalfEdgeHandle eh, std::vector<HexadHandle> &hexadVec)
	{
		HalfEdgeHH beg = handle_to_entity(eh);
		HalfEdgeHH iter = beg;
		hexadVec.push_back(beg.hedron_handle());
		while ( has_radial_half_edge(iter.handle()) && 
			   (radial_half_edge(iter.handle()).handle() != beg.mate_half_edge_handle()) )
		{
			iter = radial_half_edge(iter.handle());
			iter = mate_half_edge(iter.handle());
			hexadVec.push_back(iter.hedron_handle());
		}

		if ( !has_radial_half_edge(iter.handle()) )
		{
			iter = mate_half_edge(beg.handle());
			while( has_radial_half_edge(iter.handle()) )
			{
				iter = radial_half_edge(iter.handle());
				iter = mate_half_edge(iter.handle());
				hexadVec.push_back(iter.hedron_handle());
			}
		}
		return hexadVec.size();
	}

	/** Retrieve VertexStar in HexadMesh.
	*	\param  vh handle of the Vertex. 
	*	\param  hexadVec store the result of the VertexStar. 
	*   \return number of the Vertex in VertexStar.
	*/
	int HexadMesh::vertex_star(VertexHandle vh, std::vector<HexadHandle> &hexadVec)
	{
		Hexahedron hexad = hexahedron(vh);
		hexadVec.push_back(hexad.handle());
		for ( int i = 0; i < 6; i++ )
		{
			HalfFaceHH hf = handle_to_entity(HalfFaceHandle(hexad.first_half_face_handle().idx() + i));
			if ( hf.has_opposite_face() && has_vertex(hf.handle(), vh) )
			{
				HalfFaceHH opphf = opposite_half_face(hf.handle());
				if( std::find(hexadVec.begin(),hexadVec.end(),opphf.hedron_handle()) == hexadVec.end() )
				{
					Hexahedron oppHexad = handle_to_entity(opphf.hedron_handle());
					VertexFindIF vif(handle_to_entity(vh));
					std::vector<Vertex>::iterator iter = std::find_if(_vc.begin()+oppHexad.first_vertex_handle(), 
						                                              _vc.begin()+oppHexad.first_vertex_handle()+8,vif);
					VertexHandle newVh(iter - _vc.begin());
					vertex_star(newVh, hexadVec);
				}
			}
		}
		return hexadVec.size();
	}

	/** Retrieve PointStar in HexadMesh.
	*	\param  ph handle of the Point. 
	*	\param  hexadVec store the result of the PointStar. 
	*   \return number of the Point in PointStar.
	*/
	int HexadMesh::point_star(PointHandle ph, std::vector<HexadHandle> &hexadVec)
	{
		vertex_star(_pvc[ph.idx()], hexadVec);
		return hexadVec.size();
	}

	/** Retrieve Vertex half edge star in HexadMesh.
	*	\param  vh handle of the Vertex. 
	*	\param  halfEdgeVec store the result of the Vertex half edge Star. 
	*   \return number of the HalfEdge in halfEdgeVec.
	*/
	int HexadMesh::vertex_half_edge_star(VertexHandle vh, std::vector<HalfEdgeHandle> &halfEdgeVec)
	{
		std::vector<HexadHandle> hexVec;
		VertexFindIF vif(handle_to_entity(vh));
		std::vector<Vertex>::iterator iter;
		VertexHalfEdgeIter vhe_it;
		Vertex currV(vh);
		vertex_star(vh, hexVec);
		for (unsigned int i = 0; i < hexVec.size(); i++)
		{
			if (!is_valid(hexVec[i]))
				continue;
			iter = std::find_if(_vc.begin()+hexVec[i]*8, _vc.begin()+hexVec[i]*8+8,vif);
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
	int HexadMesh::vertex_half_face_star(VertexHandle vh, std::vector<HalfFaceHandle> &halfFaceVec)
	{
		std::vector<HexadHandle> hexVec;
		VertexFindIF vif(handle_to_entity(vh));
		std::vector<Vertex>::iterator iter;
		VertexHalfFaceIter vhf_it;
		Vertex currV(vh);
		vertex_star(vh, hexVec);
		for (unsigned int i = 0; i < hexVec.size(); i++)
		{
			if (!is_valid(hexVec[i]))
				continue;
			iter = std::find_if(_vc.begin()+hexVec[i]*8, _vc.begin()+hexVec[i]*8+8,vif);
			currV = *iter;
			for (vhf_it = vertex_half_face_iter(currV.handle()); vhf_it; ++ vhf_it)
				halfFaceVec.push_back(vhf_it.handle());
		}
		return halfFaceVec.size();
	}

	/** build opposite halfFace for HexadMesh
	*/	
	bool HexadMesh::build_opposite_half_face()
	{
		int halfFaceNum =  _npolyhedron* 6;
		int pointHandle[4];
		char buf[128];
		for ( int i=0; i<halfFaceNum; i++ )
		{
			int heIdx = i * 4;
			HalfFaceHH &hf = _hfc[i];
			pointHandle[0] = start_point_handle(HalfEdgeHandle(heIdx));
			pointHandle[1] = start_point_handle(HalfEdgeHandle(heIdx+1));
			pointHandle[2] = start_point_handle(HalfEdgeHandle(heIdx+2));
			pointHandle[3] = start_point_handle(HalfEdgeHandle(heIdx+3));
			std::sort(pointHandle, pointHandle + 4);
			sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
			if ( _ophfc.find(buf) != _ophfc.end() )
			{
				if ( _ophfc[buf]!=-1 )
				{
					HalfFaceHH &opphf = _hfc[_ophfc[buf]];
					opphf.set_opposite_face_handle(hf.handle());
					hf.set_opposite_face_handle(opphf.handle());
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
	/**
	*                       |y
	*						|
	*					    |
	*						|
	*					 010|_________________________011
	*					   /|                        /|
	*					  / |                       / |
	*					 /  |                      /  |
	*					/   |                     /   |
	*			   110 /____|____________________/111 |
	*				  |  	|                    |    |
	*				  |  000|____________________|____|________________x
	*				  |    /                     |    /001
	*				  |   /                      |   /                                  
	*				  |  /                       |  /
	*				  | /                        | /
	*              100|/_________________________|/101
	*			      /
	*				 /
	*			    /
	*			   /z
	*					 
	*             the vertex number in the HexadMesh
	*/	

	/**
	*                       |y
	*						|
	*						|
	*						|
	*					   3|_________________________2
	*					   /|                        /|
	*					  / |                       / |
	*					 /  |                      /  |
	*					/   |                     /   |
	*			     7 /____|____________________/6   |
	*				  |  	|                    |    |
	*				  |    0|____________________|____|________________x
	*				  |    /                     |    /1
	*				  |   /                      |   /                                  
	*				  |  /                       |  /
	*				  | /                        | /
	*                4|/_________________________|/5
	*			      /
	*				 /
	*			    /
	*			   /z
	*					 
	*              the vertex number of the file formate .hex
	*			    
	*			    
	*/



	/** mesh topology construct function for HexadMesh ( from vertex information )
	*/
	bool HexadMesh::build_topology()
	{
		int hfIdx, heIdx, verIdx;
		_npolyhedron = _vc.size() / 8;
		_npoint = _pc.size();
		for ( int i = 0; i < _npolyhedron; i++ )
		{
			hfIdx = i * 6;
			heIdx = i * 24;
			verIdx = i * 8;
			_hexac.push_back(Hexahedron(i));
			for ( int j = 0; j < 6; j ++ )
			{
				_hfc.push_back(HalfFaceHH(hfIdx + j));
			}
			// Halfedge of halfface 0
			_hec.push_back(HalfEdgeHH(heIdx,   verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+1, verIdx+3));
			_hec.push_back(HalfEdgeHH(heIdx+2, verIdx+7));
			_hec.push_back(HalfEdgeHH(heIdx+3, verIdx+6));
			// Halfedge of halfface 1
			_hec.push_back(HalfEdgeHH(heIdx+4, verIdx+4));
			_hec.push_back(HalfEdgeHH(heIdx+5, verIdx+5));
			_hec.push_back(HalfEdgeHH(heIdx+6, verIdx+6));
			_hec.push_back(HalfEdgeHH(heIdx+7, verIdx+7));
			// Halfedge of halfface 2
			_hec.push_back(HalfEdgeHH(heIdx+8, verIdx+1));
			_hec.push_back(HalfEdgeHH(heIdx+9, verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+10, verIdx+6));
			_hec.push_back(HalfEdgeHH(heIdx+11, verIdx+5));
			// Halfedge of halfface 3
			_hec.push_back(HalfEdgeHH(heIdx+12, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+13, verIdx+3));
			_hec.push_back(HalfEdgeHH(heIdx+14, verIdx+2));
			_hec.push_back(HalfEdgeHH(heIdx+15, verIdx+1));
			// Halfedge of halfface 4
			_hec.push_back(HalfEdgeHH(heIdx+16, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+17, verIdx+4));
			_hec.push_back(HalfEdgeHH(heIdx+18, verIdx+7));
			_hec.push_back(HalfEdgeHH(heIdx+19, verIdx+3));
			// Halfedge of halfface 5
			_hec.push_back(HalfEdgeHH(heIdx+20, verIdx+0));
			_hec.push_back(HalfEdgeHH(heIdx+21, verIdx+1));
			_hec.push_back(HalfEdgeHH(heIdx+22, verIdx+5));
			_hec.push_back(HalfEdgeHH(heIdx+23, verIdx+4));	

		}

		/**
		* build the connection between the points and the vertices
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

	/** 
	 * retrieve the Point
	 */
	Point & HexadMesh::point(PointHandle & ph)
	{
		return _pc[ph.idx()];
	}

	/** Access from vertex handle to its Point.
	*	\param  handle the index of the HalfEdgeTH, must be valid. 
	*   \return Point which is the start point of HalfEdgeTH.
	*/
	Point& HexadMesh::point(VertexHandle & handle)
	{
		return handle_to_entity(handle_to_entity(handle).point_handle());
	}


	/** Access from point iterator to its Point.
	*	\param  iter the iterator of the point, its handle must be valid. 
	*   \return Point which is the point of iterator .
	*/
	Point& HexadMesh::point(PointIter & iter)
	{		
		return handle_to_entity(iter.handle());
	}

	/** Access from vertex handle to its Point.
	*	\param  handle the index of the vertex, its handle must be valid. 
	*   \return Point which is the point of iterator.
	*/
	Point& HexadMesh::point(VertexIter & iter )
	{
		return handle_to_entity(handle_to_entity(iter.handle()).point_handle());
	}

	/** Get property of the face.
	*	\param  handle the index of the HalfFace which we want to get its property, must be valid. 
	*   \return property value.
	*/
	int HexadMesh::property(HalfFaceHH hf)
	{
		return hf.property();
	}
	
	/** Get property of the face.
	*	\param  handle the index of the HalfFace which we want to get its property, must be valid. 
	*   \return property value.
	*/
	int HexadMesh::property(HalfFaceHandle hfh)
	{
		return handle_to_entity(hfh).property();
	}
	
	/** Set property of the face.
	*	\param  handle the index of the HalfFace which we want to get its property, must be valid. 
	*/
	void HexadMesh::set_property(HalfFaceHH hf, int prop)
	{
		hf.set_property(prop);
	}

	/** Set property of the face.
	*	\param  handle the index of the HalfFace which we want to get its property, must be valid. 
	*/
	void HexadMesh::set_property(HalfFaceHandle hfh, int prop)
	{
		handle_to_entity(hfh).set_property(prop);
	}
	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve the vertex handle which the half edge is from
	*  \param handle the handle of the half edge, must be valid
	*  \return VertexHandle the vertex handle 
	*/
	VertexHandle HexadMesh::from_vertex_handle(HalfEdgeHandle hh_)
	{
		assert(hh_.is_valid());
		return handle_to_entity(hh_).start_vertex_handle();
	}

	/** Retrieve the vertex handle which the half edge is to
	*  \param handle the handle of the half edge, must be valid
	*  \return VertexHandle the vertex handle 
	*/
	VertexHandle HexadMesh::to_vertex_handle(HalfEdgeHandle hh_)
	{
		assert(hh_.is_valid());
		return next_half_edge(hh_).start_vertex_handle();
	}





	/** Retrieve Tetra stored in the TetraMesh, translate the iterator to the item
	*	\param  iter the iterator of the Tetra to retrieve, its handle must be valid. 
	*   \return Tetrahedron which indicated by the iterator
	*/
	Hexahedron & HexadMesh::iter_to_entity(HedronIter & iter)
	{
		if ( ((unsigned int)iter.handle().idx() < _hexac.size()) && (iter.handle().idx() >= 0))
		{
			return _hexac[iter.handle().idx()];
		}
		else
			return *((Hexahedron*)NULL);
	}

	/** Retrieve Vertex stored in the TetraMesh indicated by iterator
	*	\param  iter the iterator of the Vertex to retrieve, its handle must be valid. 
	*   \return Vertex which indicated by the iterator
	*/
	Vertex& HexadMesh::iter_to_entity(VertexIter& iter)
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
	Point & HexadMesh::iter_to_entity(PointIter& iter)
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
	//-----------------------------------------------------------------------------------------------------------------

	/** update the face normals
	*/
	void HexadMesh::update_face_normals()
	{
		assert(_rfn);

		_fn.resize(_hfc.size(), Vec3d(1, 0, 0));
		HFContainer::iterator it;
		HalfFaceVertexIter hfv_iter;
		Vec3d v[4];
		Vec3d n;
		unsigned int i;
		for (it = _hfc.begin(); it != _hfc.end(); ++ it)
		{
			i = 0;
			for (hfv_iter = half_face_vertex_iter((*it).handle()); hfv_iter; ++ hfv_iter)
			{
				v[i] = point(hfv_iter.handle());
				++ i;
			}

			n = (v[1] - v[0]) % (v[2] - v[0]);
			n.normalize();

			_fn[(*it).handle().idx()] = n;
		}
	}

	/** update the face normal
	*/
	Vec3d HexadMesh::update_face_normal(HalfFaceHandle & hfh_)
	{
		assert(_rfn);
		if (_fn.empty())
		{
			_fn.resize(_hfc.size(), Vec3d(1, 0, 0));
		}
		HalfFaceVertexIter hfv_iter;
		//FaceHalfedgeIter fhe_iter;
		Vec3d v[3];
		Vec3d n;
		//unsigned int i;
		//i = 0;
		hfv_iter = half_face_vertex_iter(hfh_);
		v[0] = point(hfv_iter.handle());
		++ hfv_iter;
		v[1] = point(hfv_iter.handle());
		++ hfv_iter;
		v[2] = point(hfv_iter.handle());
		//for (hfv_iter = half_face_vertex_iter(hfh_); hfv_iter; ++ hfv_iter)
		//{
		//	v[i] = point(hfv_iter.handle());
		//	++ i;
		//}
		//for (fhe_iter = face_half_edge_iter(hfh_); fhe_iter; ++ fhe_iter)
		//{
		//	v[i] = point(from_vertex_handle(fhe_iter.handle()));
		//	++ i;
		//}

		n = (v[1] - v[0]) % (v[2] - v[0]);
		n.normalize();

		_fn[hfh_.idx()] = n;
		return _fn[hfh_.idx()];
	}
	/** return the face normal
	*/
	Vec3d HexadMesh::normal(HalfFaceHandle &  hfh_)
	{
		return _fn[hfh_.idx()];
	}




	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve the first half face handle to from a tetrahedron handle
	*	\param  th_ the index of the tetrahedron, its handle must be valid. 
	*   \return HalfFaceHandle which is the handle of the half face.
	*/
	HalfFaceHandle HexadMesh::first_half_face_handle(HexadHandle th_)
	{
		return handle_to_entity(th_).first_half_face_handle();
	}

	/** Retrieve the first half edge handle to from a half handle
	*	\param  fh_ the index of the half face, its handle must be valid. 
	*   \return HalfEdgeHandle which is the handle of the half edge.
	*/
	HalfEdgeHandle HexadMesh::first_half_edge_handle(HalfFaceHandle fh_)
	{
		return handle_to_entity(fh_).first_half_edge_handle();
	}

	/** Retrieve the first half face handle to from a tetrahedron iterator
	*	\param  tIter the iterator of the tetrahedron, its handle must be valid. 
	*   \return HalfFaceHandle which is the handle of the half face.
	*/
	HalfFaceHandle HexadMesh::first_half_face_handle(HedronIter tIter_)
	{
		return iter_to_entity(tIter_).first_half_face_handle();
	}


	/** Retrieve the first vertex handle to from a half face handle
	*	\param  fh_ the handle of the half face, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle HexadMesh::first_vertex_handle(HalfFaceHandle fh_)
	{
		return handle_to_entity(handle_to_entity(fh_).first_half_edge_handle()).start_vertex_handle();
	}

	/** Retrieve the first vertex handle to from a hedron handle
	*	\param  hh_ the handle of the tetrahedron, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle HexadMesh::first_vertex_handle(PolyhedronHandle hh_)
	{
		return handle_to_entity(hh_).first_vertex_handle();
	}

	/** Retrieve the first vertex handle  from a hedron iter
	*	\param  tIter_ the iterator of the tetrahedron, its handle must be valid. 
	*   \return VertexHandle which is the handle of vertex.
	*/
	VertexHandle HexadMesh::first_vertex_handle(HedronIter tIter_)
	{
		return handle_to_entity(tIter_.handle()).first_vertex_handle();
	}

	/** Retrieve the next vertex handle from a hedron handle 
	*  \param hh_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle HexadMesh::next_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && hh_.is_valid());
		return VertexHandle ((cvh_.idx() + 1) % 8 + handle_to_entity(hh_).first_vertex_handle().idx()); //((int)(cvh_.idx() / 4)) * 4
	}

	/** Retrieve the next vertex handle from a hedron iterator 
	*  \param tIter_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle HexadMesh::next_vertex_handle(HedronIter tIter_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && tIter_.handle().is_valid());
		return VertexHandle((tIter_.handle().idx() + 1) % 8 + handle_to_entity(tIter_.handle()).first_vertex_handle().idx());//((int)(hh_.idx() / 4)) * 4
	}

	/** Retrieve the prev vertex handle from a hedron handle 
	*  \param hh_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle HexadMesh::prev_vertex_handle(PolyhedronHandle hh_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && hh_.is_valid());
		return VertexHandle ((cvh_.idx() - 1) % 8 + handle_to_entity(hh_).first_vertex_handle().idx()); //((int)(cvh_.idx() / 4)) * 4
	}

	/** Retrieve the prev vertex handle from a hedron iterator 
	*  \param tIter_ the polyhedral handle
	*  \param cvh_ the current vertex handle
	*  \return VertexHandle which is the next vertex handle
	*/
	VertexHandle HexadMesh::prev_vertex_handle(HedronIter tIter_, VertexHandle cvh_)
	{
		assert(cvh_.is_valid() && tIter_.handle().is_valid());
		return VertexHandle((tIter_.handle().idx() - 1) % 8 + handle_to_entity(tIter_.handle()).first_vertex_handle().idx());//((int)(hh_.idx() / 4)) * 4
	}

	/** Check whether the HalfFaceTH is boundary half face.
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::is_boundary(HalfFaceHandle & hfh_)
	{
		return !handle_to_entity(hfh_).has_opposite_face();
	}

	/** Check whether the vertex handle  is boundary vertex.
	*	\param  vh_ the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::is_boundary(VertexHandle & vh_)
	{
		return _brv[handle_to_entity(vh_).point_handle().idx()];
		//return !(_bv.find(handle_to_entity(vh_).point_handle()) == _bv.end());
	}
	/** Check whether the hedron is boundary hedron, if one of its vertices is boundary point, it is a boundary hedron
	*	\param  handle the index of the HalfFaceTH, must be valid. 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::is_boundary(PolyhedronHandle & hh_)
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
		//return is_boundary()
	}

	/** Check whether the hedron is boundary hedron
	*	\param  handle the index of the HalfFaceTH, must be valid.
	*   \param type_ the boundary type. 0x0001 for half face boundary while 0x0002 for only vertex boundary 
	*   \return true  yes.
	*           false no
	*/
	bool HexadMesh::is_boundary(HedronHandle & hh_, unsigned int type_)
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

	//-----------------------------------------------------------------------------------------------------------------

	/** Build the boundary vertices of the mesh
	*/
	void HexadMesh::update_boundary_point()
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


	// clean _free_hedron_list and _free_point_list & filling the vacant
	void HexadMesh::clean_garbage()
	{
		// clean _free_tetra_list
		for (unsigned int i = 0; i < _free_hedron_list.size(); i++)
		{
			if (_free_hedron_list[i].idx() != _npolyhedron - 1)
			{
				Hexahedron lastHexad = handle_to_entity(HexadHandle(_npolyhedron - 1));
				HalfFaceHH lastHf = handle_to_entity(lastHexad.first_half_face_handle());
				HalfFaceHH hf = handle_to_entity(handle_to_entity(_free_hedron_list[i]).first_half_face_handle());

				//reset last halfface's opposite relation to new position
				for ( int j = 0; j < 6; j ++ )
				{
					if ( lastHf.has_opposite_face() )
					{
						HalfFaceHH& opphf = handle_to_entity(lastHf.opposite_face_handle());
						opphf.set_opposite_face_handle(HalfFaceHandle(_free_hedron_list[i].idx() * 6 + j));
						_hfc[hf.handle()].set_opposite_face_handle(opphf.handle());
					}
					else
					{
						_hfc[hf.handle()].set_opposite_face_handle(HalfFaceHandle(-1));

						// reset _ophfc' Map relation for the last Hexad
						char buf[128];
						PointHandle pointHandle[4];
						HalfEdgeHandle lastHe = lastHf.first_half_edge_handle();
						pointHandle[0] = start_point_handle(HalfEdgeHandle(lastHe));
						pointHandle[1] = start_point_handle(HalfEdgeHandle(lastHe+1));
						pointHandle[2] = start_point_handle(HalfEdgeHandle(lastHe+2));
						pointHandle[3] = start_point_handle(HalfEdgeHandle(lastHe+3));
						std::sort(pointHandle, pointHandle+4);
						sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
						_ophfc[buf] = _free_hedron_list[i]*6+j;
					}
					hf = next_half_face(hf);
					lastHf = next_half_face(lastHf);
				}

				//reset point handle 
				for ( int j=0; j<8; j++ )
				{
					Vertex lastV = handle_to_entity(VertexHandle(lastHexad.first_vertex_handle().idx()+j));
					Vertex &v = handle_to_entity(VertexHandle(handle_to_entity(_free_hedron_list[i]).first_vertex_handle().idx()+j));
					v.set_point_handle(lastV.point_handle());
				}
			}

			//delete all the topology relation
			_hexac.pop_back();
			for ( int j=0; j<6; j++ )
				_hfc.pop_back();
			for ( int j=0;j<24;j++ )
				_hec.pop_back();
			for ( int j=0; j<8; j++ )
				_vc.pop_back();
			_npolyhedron = _hexac.size();
		}

		// clean _free_point_list
		for (unsigned int i = 0; i < _free_point_list.size(); i++)
		{
			// This point has no reference, should be delete. 2 rule.
			// #1 We exchange it with the last point in vector, then pop vector.
			// #2 The erase point is the last point in vector, just pop it.
			// #3 Revise the opposite mapping relation
			if ( _free_point_list[i].idx() != _npoint-1 )
			{
				for (unsigned int j=0; j < _hfc.size(); j ++)
				{
					if ( has_point(HalfFaceHandle(j), PointHandle(_npoint-1)) )
					{
						char buf[128];
						int tmpHf;

						//Revise opposite mapping relation
						PointHandle pointHandle[4];
						pointHandle[0] = start_point_handle(HalfEdgeHandle(j*4));
						pointHandle[1] = start_point_handle(HalfEdgeHandle(j*4+1));
						pointHandle[2] = start_point_handle(HalfEdgeHandle(j*4+2));
						pointHandle[3] = start_point_handle(HalfEdgeHandle(j*4+3));
						std::sort(pointHandle, pointHandle+4);
						sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
						if ( _ophfc.find(buf)!=_ophfc.end() )
						{
							tmpHf = _ophfc[buf];
							_ophfc.erase(buf);

							pointHandle[3] = _free_point_list[i];
							std::sort(pointHandle, pointHandle+4);
							sprintf_s(buf,"%d %d %d %d",pointHandle[0],pointHandle[1],pointHandle[2],pointHandle[3]);
							_ophfc[buf] = tmpHf;
						}
					}
				}

				//Revise PointHandle in Vertex
				std::vector<Vertex>::iterator iter;
				VertexFindIF vif = Vertex(_npoint-1);
				while ( (iter=find_if(_vc.begin(),_vc.end(),vif))!=_vc.end() )
				{
					(*iter).set_point_handle(_free_point_list[i]);
				}
			}

			//--------------------------new added---------------------------------------//
			_pvc[_free_point_list[i].idx()] = _pvc[_npoint - 1];
			_pvc.pop_back();
			//--------------------------end added---------------------------------------//
			_pc[_free_point_list[i]] = _pc[_npoint-1];
			_pc.pop_back();
			_npoint--;
		}
	}



//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------
