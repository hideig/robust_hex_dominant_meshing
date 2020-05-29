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
/**********************************************************************************************************************
*                                                                                                                     *
*                               VolumeMesh                                                                            *
*      Copyright (C) 2010 by Computer Aided Designed Group                                                            *
*        State Key Lab. of CAD & CG, Zhejiang University                                                              *
*         This project is created by Chuhua Xian, 2010.2                                                              *
*                     Email: chuhuaxian@gmail.com																	  *	
*         Modified by Xiaoshen Chen, 2010.03																		  *
*					   Email: chinimei@163.com																		  *
***********************************************************************************************************************/ 
/*---------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		       *
*					   Email: chuhuaxian@gmail.com																       *
*---------------------------------------------------------------------------------------------------------------------*/ 



#ifndef _VOLUME_MESH_TOPOLOGY_H_
#define _VOLUME_MESH_TOPOLOGY_H_

//---------------------------------------------------------------------------------------------------------------------

#include <VolumeMesh/Mesh/Handles.h>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

	typedef PolyhedronHandle TetraHandle;
	typedef PolyhedronHandle HexadHandle;

	/** 
	 * definition of Polyhedron, it is the basic class
	 */
	class Polyhedron
	{
	public:
		Polyhedron()
		{
		}

		Polyhedron(unsigned int type) : _pt(type)
		{

		}

		virtual ~Polyhedron()
		{
		}
	public:
		unsigned int polyhedron_type()
		{
			return _pt;
		}
	protected:
		unsigned int _pt; // 0x0001 for tetrahedron, 0x0002 for hexahedron
	};


//---------------------------------------------------------------------------------------------------------------------
	
	/// definition of Tetrahedron
	class Tetrahedron : public Polyhedron
	{
		
	public:

		Tetrahedron() : Polyhedron(0x0001)
		{
		};
		Tetrahedron(int idx) : _tetra(idx)
		{
		};
		~Tetrahedron()
		{
		}

	public:
		TetraHandle handle();
		HalfFaceHandle first_half_face_handle();
		VertexHandle first_vertex_handle();
	protected:
		TetraHandle _tetra;
	};

	/// definition of half face
	class HalfFace
	{
	public:
		HalfFace()
		{
			_ft = 0;
		}

		HalfFace(unsigned int type) : _ft (type)
		{
		}

		virtual ~ HalfFace()
		{
		}
	public:
		unsigned int face_type()
		{
			return _ft;
		}
	protected:
		unsigned int _ft; /**< the type of the face, 0x0001 for triangle face, 0x0002 for quad face */
	};

	/// definition of halfedge
	class HalfEdge
	{
	public:
		HalfEdge()
		{
			_et = 0;
		}

		HalfEdge(unsigned int type) : _et (type)
		{
		}

		virtual ~ HalfEdge()
		{
		}
	public:
		unsigned int edge_type()
		{
			return _et;
		}
	protected:
		unsigned int _et;
	};


	/// definition of HalfFaceTH
	class HalfFaceTH : public HalfFace
	{
	public:
		HalfFaceTH(int hf,int opphf = -1) : _hf(hf), _opphf(opphf), HalfFace(0x0001)
		{
			_prop = -1;
		};

		HalfFaceTH(HalfFaceHandle &phf, HalfFaceHandle &popphf = HalfFaceHandle(-1)) : _hf(phf), _opphf(popphf), HalfFace(0x0001)
		{
			_prop = -1;
		};

		~HalfFaceTH()
		{
		}
	public:
		HalfFaceHandle handle();
		HalfFaceHandle opposite_face_handle();
		bool has_opposite_face();
		bool set_opposite_face_handle(HalfFaceHandle hf);

		HalfFaceHandle next_half_face_handle();
		HalfFaceHandle prev_half_face_handle();
		HalfFaceHandle mid_half_face_handle();

		const HalfFaceHandle next_half_face_handle() const;
		const HalfFaceHandle prev_half_face_handle() const;
		const HalfFaceHandle mid_half_face_handle()  const;

		HalfEdgeHandle first_half_edge_handle();
		TetraHandle hedron_handle();

		int property();
		void set_property(int prop);
	private:
		HalfFaceHandle _hf;
		HalfFaceHandle _opphf;	//opposite half face handle
		int _prop;        /**< the property of the face */
	};

	/// definition of HalfEdgeTH
	class HalfEdgeTH : public HalfEdge
	{
	public:
		HalfEdgeTH(int iHE, int iV) : _he(iHE), _vs(iV), HalfEdge(0x0001)
		{
		};
		HalfEdgeTH(HalfEdgeHandle & phe, VertexHandle & pv) : _he(phe),_vs(pv), HalfEdge(0x0001)
		{
		};
		~ HalfEdgeTH()
		{
		}

	public:
		HalfEdgeHandle handle();
		TetraHandle hedron_handle();
		HalfFaceHandle half_face_handle();
		HalfEdgeHandle mate_half_edge_handle();
		HalfEdgeHandle next_half_edge_handle();
		HalfEdgeHandle prev_half_edge_handle();
		const HalfEdgeHandle next_half_edge_handle() const;
		const HalfEdgeHandle prev_half_edge_handle() const;
		VertexHandle start_vertex_handle();
	private:
		HalfEdgeHandle _he;
		VertexHandle _vs;		//global
	};

	/// function object to find Tetrahedron
	class TetrahedronFindIF
	{
	public:
		TetrahedronFindIF(Tetrahedron t):_t(t)
		{

		};
		~ TetrahedronFindIF()
		{
		}
	public:
		bool operator ()(Tetrahedron& t);
		void operator = (Tetrahedron& t);
	private:
		Tetrahedron _t;
	};

	/// function object to find TetraHandle
	class TetraHandleFindIF
	{
	public:
		TetraHandleFindIF(TetraHandle th) : _th(th)
		{
		};
	public:
		bool operator ()(TetraHandle& th);
	private:
		TetraHandle _th;
	};

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------



/// namespace Hexad
//namespace Hexad
//{
//---------------------------------------------------------------------------------------------------------------------

	/// definition of hexahedron
	class Hexahedron : public Polyhedron
	{
	public:
		Hexahedron() : Polyhedron(0x0002)
		{
		};
		Hexahedron(int idx) : _hexad(idx), Polyhedron(0x0002)
		{
		};
		~Hexahedron()
		{
		}
	public:
		HexadHandle handle();
		HalfFaceHandle first_half_face_handle();
		VertexHandle first_vertex_handle();
	protected:
	private:
		HexadHandle _hexad;
	};

	/// definition of HalfFaceTH
	class HalfFaceHH : public HalfFace
	{
	public:
		HalfFaceHH() : HalfFace(0x0002)
		{
			_prop = -1;
		}
		HalfFaceHH(int hf,int opphf = -1): _hf(hf), _opphf(opphf), HalfFace(0x0002)
		{
			_prop = -1;
		};
		HalfFaceHH(HalfFaceHandle &phf, HalfFaceHandle &popphf = HalfFaceHandle(-1)) : _hf(phf), _opphf(popphf), HalfFace(0x0002)
		{
			_prop = -1;
		};
		~HalfFaceHH()
		{
		}

	public:
		HalfFaceHandle handle();
		HalfFaceHandle opposite_face_handle();
		HalfFaceHandle opposite_face_handle() const;
		bool has_opposite_face();
		bool set_opposite_face_handle(HalfFaceHandle hf);
		HalfFaceHandle next_half_face_handle();
		HalfFaceHandle prev_half_face_handle();
		HalfFaceHandle next_half_face_handle() const;
		HalfFaceHandle prev_half_face_handle() const;
		HalfEdgeHandle first_half_edge_handle();
		HexadHandle hedron_handle();
		int property();
		void set_property(int prop);
	protected:
	private:
		HalfFaceHandle _hf;
		HalfFaceHandle _opphf;	//opposite HalfFaceTH handle
		int _prop;              /**< the property of the face */
	};

	/// definition of HalfEdgeTH
	class HalfEdgeHH : public HalfEdge
	{
	public:
		HalfEdgeHH() : HalfEdge(0x0002)
		{

		}
		HalfEdgeHH(int iHE, int iV) : _he(iHE), _vs(iV), HalfEdge(0x0002)
		{
		};
		HalfEdgeHH(HalfEdgeHandle &phe, VertexHandle &pv) : _he(phe), _vs(pv), HalfEdge(0x0002)
		{
		};
	public:
		HexadHandle hedron_handle();
		HalfFaceHandle half_face_handle();
		HalfEdgeHandle mate_half_edge_handle();
		HalfEdgeHandle mate_half_edge_handle() const;
		HalfEdgeHandle handle();
		HalfEdgeHandle next_half_edge_handle();
		HalfEdgeHandle prev_half_edge_handle();
		HalfEdgeHandle next_half_edge_handle() const;
		HalfEdgeHandle prev_half_edge_handle() const;
		VertexHandle start_vertex_handle();
	protected:
	private:
		HalfEdgeHandle _he;
		VertexHandle _vs;		//global
	};

	/// function object to find Hexahedron
	class HexahedronFindIF
	{
	public:
		HexahedronFindIF(Hexahedron h) : _h(h)
		{
		};
		~HexahedronFindIF()
		{

		}
	public:
		//bool operator ()(Hexahedron& h);
		//void operator = (Hexahedron& h);
		/** Function object to find Hexahedron.
		*   \return true  found.
		*			false not found.
		*/
		bool operator () (Hexahedron& h)
		{
			return _h.handle() == h.handle();
		}

		/** Assignment function of Function object to find Hexahedron.
		*/
		void operator = (Hexahedron& h)
		{
			_h = h;
		}

	

	private:
		Hexahedron _h;
	};

	/// function object to find HexadHandle
	class HexadHandleFindIF
	{
	public:
		HexadHandleFindIF(HexadHandle hh):_hh(hh)
		{
		};
		~HexadHandleFindIF()
		{
		}
	public:
		/** Function object to find HexaHandle.
		*   \return true  found.
		*			false not found.
		*/
		bool operator ()(HexadHandle& hh)
		{
			return _hh == hh;
		}
	private:
		HexadHandle _hh;
	};

//---------------------------------------------------------------------------------------------------------------------
/*}	// namespace Hexad
//---------------------------------------------------------------------------------------------------------------------
*/	
	/// definition of Vertex
	class Vertex
	{
	public:
		Vertex(PointHandle &pp) : _p(pp), _next(-1)
		{
		};
		Vertex(int iP) : _p(iP), _next(-1)
		{
		};
		Vertex(VertexHandle &vh, PointHandle &pp) : _vh(vh), _p(pp), _next(-1)
		{
		};
		Vertex(int vh, int pp):_vh(vh), _p(pp), _next(-1)
		{
		};
	public:
		VertexHandle handle()
		{
			return _vh;
		}
		VertexHandle next_vertex_handle()
		{
			return _next;
		}
		void set_next_vertex_handle(VertexHandle & handle)
		{
			_next = handle;
		}
		PointHandle point_handle()
		{
			return _p; 
		};
		void set_point_handle(PointHandle& handle) 
		{
			_p = handle;
		}
	private:
		VertexHandle _vh;
		VertexHandle _next;
		PointHandle _p;
		//bool _is_boundary; /**< if the vertex is boundary vertex */
	};	

	/// function object to find Vertex
	class VertexFindIF
	{
	public:
		VertexFindIF(Vertex v) : _v(v)
		{
		};
		// bool operator ()(Vertex& v);
		// void operator = (Vertex& v);
	
		/** Function object to find equal Vertex.
		*   \return true  found.
		*			false not found.
		*/
		bool operator ()(VolumeMesh::Vertex &v)
		{
			return _v.point_handle()==v.point_handle();
		}

		/** Assignment function of Function object to find Vertex.
		*/
		void operator =(VolumeMesh::Vertex &v)
		{
			_v = v;
		}

	private:
		Vertex _v;
	};

	/// function object to find VertexHandle
	class VertexHandleFindIF
	{
	public:
		VertexHandleFindIF(VertexHandle vh):_vh(vh)
		{

		};
		//bool operator ()(VertexHandle& vh);

		/** Function object to find equal Vertex.
		*   \return true  found.
		*			false not found.
		*/
		bool operator ()(VolumeMesh::VertexHandle& vh)
		{
			return _vh == vh;
		}
	private:
		VertexHandle _vh;
	};

//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------


#endif