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


#include <VolumeMesh/Mesh/Topology.h>

namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------

	extern const int TNEXT[4][4] = 
	{
		-1,  2,  3,  1,
		 3, -1,  0,  2,
		 1,  3, -1,  0,
		 2,  0,  1,  -1
	};

	extern const int TPREV[4][4] = 
	{
		-1,  3,  1,  2,
		 2, -1,  3,  0,
		 3,  0, -1,  1,
		 1,  2,  0, -1
	};

	extern const int THEPosition[4][4] = 
	{
		-1,  0,  1,  2,
		 0, -1,  1,  2,
		 0,  1, -1,  2,
		 0,  1,  2,  -1
	};

	/** Retrieve TetraHandle.
	*   \return TetraHandle.
	*/
	TetraHandle Tetrahedron::handle()
	{
		return _tetra;
	}

	/** Retrieve Tetrahedron's first HalfFaceTH handle.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle Tetrahedron::first_half_face_handle()
	{
		return HalfFaceHandle(_tetra.idx() << 2); 
	}

	/** Retrieve Tetrahedron's first Vertex.
	*   \return VertexHandle.
	*/
	VertexHandle Tetrahedron::first_vertex_handle()
	{
		return VertexHandle(_tetra.idx()<<2);
	}
	
	/** Retrieve HalfFaceHandle.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceTH::handle()
	{
		return _hf;
	}

	/** Retrieve Opposite HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceTH::opposite_face_handle()
	{
		return _opphf;
	}

	/** Check whether HalfFaceTH has opposite HalfFaceTH.
	*   \return true  yes
	*			false no
	*/
	bool HalfFaceTH::has_opposite_face()
	{
		return _opphf.idx() != -1;
	}

	/** Set HalfFaceTH's Opposite HalfFaceTH relation.
	*	\param hf opposite HalfFaceTH handle
	*   \return true success.
	*			false fail.
	*/
	bool HalfFaceTH::set_opposite_face_handle(HalfFaceHandle hf)
	{
		this->_opphf = hf;
		return true;
	}

	/** Retrieve next HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceTH::next_half_face_handle()
	{
		return HalfFaceHandle(((_hf.idx() >> 2) << 2) + (_hf.idx() + 1) % 4);
	}

	/** Retrieve previous HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceTH::prev_half_face_handle()
	{
		return HalfFaceHandle(((_hf.idx() >> 2) << 2) + (_hf.idx() + 3) % 4);
	}

	/** Retrieve middle HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceTH::mid_half_face_handle()
	{
		return HalfFaceHandle(
			((_hf.idx() >> 2) << 2) + (_hf.idx() + 2) % 4
			);
	}

	/** Retrieve next HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	const HalfFaceHandle HalfFaceTH::next_half_face_handle() const
	{
		return HalfFaceHandle(((_hf.idx() >> 2) << 2) + (_hf.idx() + 1) % 4);
	}

	/** Retrieve previous HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	const HalfFaceHandle HalfFaceTH::prev_half_face_handle() const
	{
		return HalfFaceHandle(((_hf.idx() >> 2) << 2) + (_hf.idx() + 3) % 4);
	}

	/** Retrieve middle HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	const HalfFaceHandle HalfFaceTH::mid_half_face_handle() const
	{
		return HalfFaceHandle(
			((_hf.idx() >> 2) << 2) + (_hf.idx() + 2) % 4
			);
	}

	/** Retrieve HalfFaceTH's first HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfFaceTH::first_half_edge_handle()
	{
		return HalfEdgeHandle(_hf.idx() * 3);
	}

	/** Retrieve TetraHandle which HalfFaceTH belongs to.
	*   \return TetraHandle.
	*/
	TetraHandle HalfFaceTH::hedron_handle()
	{
		return TetraHandle(_hf.idx()>>2);
	}

	int HalfFaceTH::property()
	{
		return _prop;
	}

	void HalfFaceTH::set_property(int prop)
	{
		_prop = prop;
	}

	//-----------------------------------------------------------------------------------------------------------------

	/** Retrieve TetraHandle which HalfEdgeTH belongs to.
	*   \return TetraHandle.
	*/
	TetraHandle HalfEdgeTH::hedron_handle()
	{
		return TetraHandle(_he.idx()/12);
	}

	/** Retrieve HalfFaceHandle which HalfEdgeTH belongs to.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfEdgeTH::half_face_handle()
	{
		return HalfFaceHandle(_he.idx() / 3);
	}

	/** Retrieve HalfEdgeHandle.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeTH::handle()
	{
		return this->_he;
	}

	/** Retrieve HalfEdgeTH's Mate HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeTH::mate_half_edge_handle()
	{
		int mateHf = TPREV[(_he.idx() / 3) % 4][_vs.idx() % 4];
		int mateV  = TNEXT[(_he.idx() / 3) % 4][_vs.idx() % 4];
		return HalfEdgeHandle(
			((_he.idx()/12) * 4 + mateHf) * 3 + THEPosition[mateHf][mateV]
			);
	}

	/** Retrieve HalfEdgeTH's next HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeTH::next_half_edge_handle()
	{
		int nextVertexOfHE =  TNEXT[(_he.idx() / 3) % 4][_vs.idx() % 4];
		return HalfEdgeHandle(
			                 (_he.idx() / 3) * 3 + THEPosition[(_he.idx() / 3) % 4][nextVertexOfHE]
		                     );
	}

	/** Retrieve HalfEdgeTH's previous HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeTH::prev_half_edge_handle()
	{
		int prevVertexOfHE = TPREV[(_he.idx() / 3) % 4][_vs.idx() % 4];
		return HalfEdgeHandle(
			(_he.idx() / 3) * 3 + THEPosition[(_he.idx() / 3) % 4][prevVertexOfHE]
			);
	}

	/** Retrieve HalfEdgeTH's next HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	const HalfEdgeHandle HalfEdgeTH::next_half_edge_handle() const
	{
		int nextVertexOfHE =  TNEXT[(_he.idx() / 3) % 4][_vs.idx() % 4];
		return HalfEdgeHandle(
			(_he.idx() / 3) * 3 + THEPosition[(_he.idx() / 3) % 4][nextVertexOfHE]
		);
	}

	/** Retrieve HalfEdgeTH's previous HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	const HalfEdgeHandle HalfEdgeTH::prev_half_edge_handle() const
	{
		int prevVertexOfHE = TPREV[(_he.idx() / 3) % 4][_vs.idx() % 4];
		return HalfEdgeHandle(
			(_he.idx() / 3) * 3 + THEPosition[(_he.idx() / 3) % 4][prevVertexOfHE]
		);
	}

	/** Retrieve HalfEdgeTH's start Vertex.
	*   \return VertexHandle.
	*/
	VertexHandle HalfEdgeTH::start_vertex_handle()
	{
		return _vs;
	}

	/** Function object to find TetraHedron.
	*   \return true  found.
	*			false not found.
	*/
	bool TetrahedronFindIF::operator ()(VolumeMesh::Tetrahedron &t)
	{
		return _t.handle() == t.handle();
	}

	/** Assignment function of Function object to find TetraHedron.
	*/
	void TetrahedronFindIF::operator =(VolumeMesh::Tetrahedron &t)
	{
		_t = t;
	}

	/** Function object to find TetraHedron.
	*   \return true  found.
	*			false not found.
	*/
	bool TetraHandleFindIF::operator ()(VolumeMesh::TetraHandle& th)
	{
		return _th == th;
	}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------

//namespace Hexad
//{
////---------------------------------------------------------------------------------------------------------------------
//	
	

	/**
	*                       |y
	*						|
	*						|
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
	*			    
	*					 
    *
	*/

							


	
	extern const int HMate[6][8] = 
	{
		-1, -1,  3,  4, -1, -1,  2,  1,
		-1, -1, -1, -1,  5,  2,  0,  4,
		-1,  3,  0, -1, -1,  5,  1, -1,
		 4,  5,  2,  0, -1, -1, -1, -1,
		 5, -1, -1,  3,  1, -1, -1,  0,
		 3,  2, -1, -1,  4,  1, -1, -1
	};

	extern const int HHEPosition[6][8] = 
	{
		-1, -1,  0,  1, -1, -1,  3,  2,
		-1, -1, -1, -1,  0,  1,  2,  3,
		-1,  0,  1, -1, -1,  3,  2, -1,
		 0,  3,  2,  1, -1, -1, -1, -1,
		 0, -1, -1,  3,  1, -1, -1,  2,
		 0,  1, -1, -1,  3,  2, -1, -1
	};

	extern const int HNEXT[6][8] =
	{
		-1, -1,  3,  7, -1, -1,  2,  6,
		-1, -1, -1, -1,  5,  6,  7,  4,
		-1,  2,  6, -1, -1,  1,  5, -1,
		 3,  0,  1,  2, -1, -1, -1, -1,
		 4, -1, -1,  0,  7, -1, -1,  3,
		 1,  5, -1, -1,  0,  4, -1, -1
	};

	/** Retrieve HexadHandle.
	*   \return HexadHandle.
	*/
	HexadHandle Hexahedron::handle()
	{
		return _hexad;
	}

	/** Retrieve Hexahedron's first HalfFaceTH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle Hexahedron::first_half_face_handle()
	{
		return HalfFaceHandle(_hexad * 6);
	}

	/** Retrieve Hexahedron's first Vertex.
	*   \return VertexHandle.
	*/
	VertexHandle Hexahedron::first_vertex_handle()
	{
		return VertexHandle(_hexad << 3);
	}

	/** Retrieve HalfFaceHandle.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::handle()
	{
		return _hf;
	}

	/** Retrieve Opposite HalfFaceHH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::opposite_face_handle()
	{
		return _opphf;
	}

	/** Retrieve Opposite HalfFaceHH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::opposite_face_handle() const
	{
		return _opphf;
	}

	/** Check whether HalfFaceHH has opposite HalfFaceHH.
	*   \return true  yes
	*			false no
	*/
	bool HalfFaceHH::has_opposite_face()
	{
		return _opphf.idx()!=-1;
	}

	/** Set HalfFaceHH's Opposite HalfFaceHH relation.
	*	\param hf opposite HalfFaceTH handle
	*   \return true success.
	*			false fail.
	*/
	bool HalfFaceHH::set_opposite_face_handle(HalfFaceHandle hf)
	{
		_opphf = hf;
		return true;
	}

	/** Retrieve next HalfFaceHH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::next_half_face_handle()
	{
		return HalfFaceHandle(_hf / 6 * 6 + (_hf + 1) % 6);
	}

	/** Retrieve previous 
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::prev_half_face_handle()
	{
		return HalfFaceHandle(_hf / 6 * 6 + (_hf + 5) % 6);
	}

	/** Retrieve next HalfFaceHH.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::next_half_face_handle() const
	{
		return HalfFaceHandle(_hf / 6 * 6 + (_hf + 1) % 6);
	}

	/** Retrieve previous 
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfFaceHH::prev_half_face_handle() const
	{
		return HalfFaceHandle(_hf / 6 * 6 + (_hf + 5) % 6);
	}

	/** Retrieve HalfFaceHH's first HalfEdgeHH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfFaceHH::first_half_edge_handle()
	{
		return HalfEdgeHandle(_hf << 2);
	}

	/** Retrieve HexadHandle which HalfFaceHH belongs to.
	*   \return HexadHandle.
	*/
	HexadHandle HalfFaceHH::hedron_handle()
	{
		return HexadHandle(_hf / 6);
	}

	int HalfFaceHH::property()
	{
		return _prop;
	}

	void HalfFaceHH::set_property(int prop)
	{
		_prop = prop;
	}

	/** Retrieve HexadHandle which HalfEdgeHH belongs to.
	*   \return HexadHandle.
	*/
	HexadHandle HalfEdgeHH::hedron_handle()
	{
		return HexadHandle(_he / 24);
	}

	/** Retrieve HalfFaceHandle which HalfEdgeHH belongs to.
	*   \return HalfFaceHandle.
	*/
	HalfFaceHandle HalfEdgeHH::half_face_handle()
	{
		return HalfFaceHandle(_he >> 2);
	}

	/** Retrieve HalfEdgeHH's start Vertex.
	*   \return VertexHandle.
	*/
	VertexHandle HalfEdgeHH::start_vertex_handle()
	{
		return _vs;
	}

	/** Retrieve HalfEdgeHH's Mate HalfEdgeHH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::mate_half_edge_handle()
	{
		int mateHF = HMate[half_face_handle() % 6][start_vertex_handle() % 8];
		int vertex = HNEXT[half_face_handle() % 6][start_vertex_handle() % 8];
		int hePos  = HHEPosition[mateHF][vertex];
		return HalfEdgeHandle( _he / 24 * 24 + mateHF * 4 + hePos );
	}

	///** Retrieve HalfEdgeHH's Mate HalfEdgeHH.
	//*   \return HalfEdgeHandle.
	//*/
	//HalfEdgeHandle HalfEdgeHH::mate_half_edge_handle() const
	//{
	//	int mateHF = HMate[half_face_handle() % 6][start_vertex_handle() % 8];
	//	int vertex = HNEXT[half_face_handle() % 6][start_vertex_handle() % 8];
	//	int hePos  = HHEPosition[mateHF][vertex];
	//	return HalfEdgeHandle( _he / 24 * 24 + mateHF * 4 + hePos );
	//}


	/** Retrieve HalfEdgeHandle.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::handle()
	{
		return _he;
	}

	/** Retrieve HalfEdgeHH's next HalfEdgeHH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::next_half_edge_handle()
	{
		return HalfEdgeHandle(_he / 4 * 4 + (_he + 1) % 4);
	}

	/** Retrieve HalfEdgeHH's previous HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::prev_half_edge_handle()
	{
		return HalfEdgeHandle(_he/4 * 4 + (_he + 3)%4);
	}

	/** Retrieve HalfEdgeHH's next HalfEdgeHH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::next_half_edge_handle() const
	{
		return HalfEdgeHandle(_he / 4 * 4 + (_he + 1) % 4);
	}

	/** Retrieve HalfEdgeHH's previous HalfEdgeTH.
	*   \return HalfEdgeHandle.
	*/
	HalfEdgeHandle HalfEdgeHH::prev_half_edge_handle() const
	{
		return HalfEdgeHandle(_he/4 * 4 + (_he + 3)%4);
	}


////---------------------------------------------------------------------------------------------------------------------
//}	// namespace Hexad
//---------------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------
