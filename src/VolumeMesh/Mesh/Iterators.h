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
* this file defines the iterators for vertex, point, polyhedron                                                       *
* $ written by xianc $                                                                                                *
* $ Version 0.1 $                                                                                                     *
* $ Date: 2010-02-08 $                                                                                                *
**********************************************************************************************************************/
/*---------------------------------------------------------------------------------------------------------------------*
*         Modified by Chuhua Xian, 2010.06																		       *
*					   Email: chuhuaxian@gmail.com																       *
*---------------------------------------------------------------------------------------------------------------------*/ 

#ifndef _VOLUME_MESH_ITERATOR_
#define _VOLUME_MESH_ITERATOR_


#include <VolumeMesh/Mesh/Handles.h>



namespace VolumeMesh
{

	/** A VertexIterator is used to iterate over all vertex in the mesh,
	from begin() to end(), just like fr STL containers.  Its
	operator* returns a index.
	*/
	class VertexIterator
	{
	public:
		/** Constructor. 
		\internal
		*/
		VertexIterator() : idx_(-1)
		{
		}

		VertexIterator(unsigned int _idx) : idx_(_idx) 
		{
		}
		~VertexIterator()
		{
		}
		/// Get  index from iterator
		unsigned int & operator*()  
		{ 
			return idx_; 
		}
		/// Get index from iterator
		unsigned int* operator->() 
		{
			return &idx_; 
		}
		/// Comparison
		bool operator==(const VertexIterator & _rhs) const 
		{ 
			return idx_==_rhs.idx_;
		}
		/// Comparison
		bool operator!=(const VertexIterator& _rhs) const 
		{
			return !(*this == _rhs); 
		}
		/// Pre-increment
		VertexIterator & operator ++ () 
		{ 
			++idx_; 
			return *this; 
		}
		/// Pre-increment
		VertexIterator & operator -- () 
		{ 
			-- idx_; 
			return *this; 
		}

		/// handle
		VertexHandle handle()
		{
			return VertexHandle(idx_);
		}
	private:
		unsigned int idx_;
	};

	/** A PointIterator is used to iterate over all vertex in the mesh,
	from begin() to end(), just like fr STL containers.  Its
	operator* returns a index.
	*/
	class PointIterator
	{
	public:
		/** Constructor. 
		\internal
		*/
		PointIterator() : idx_(-1)
		{
		}
		PointIterator(unsigned int _idx) : idx_(_idx) 
		{
		}
		~PointIterator()
		{
		}
		/// Get  index from iterator
		unsigned int & operator*()  
		{ 
			return idx_; 
		}
		/// Get index from iterator
		unsigned int* operator->() 
		{
			return &idx_; 
		}
		/// Comparison
		bool operator==(const PointIterator & _rhs) const 
		{ 
			return idx_==_rhs.idx_;
		}
		/// Comparison
		bool operator!=(const PointIterator& _rhs) const 
		{
			return !(*this == _rhs); 
		}
		/// Pre-increment
		PointIterator & operator ++ () 
		{ 
			++idx_; 
			return *this; 
		}
		/// Pre-increment
		PointIterator & operator -- () 
		{ 
			-- idx_; 
			return *this; 
		}

		/// handle
		PointHandle handle()
		{
			return PointHandle(idx_);
		}
	private:
		unsigned int idx_;
	};



	/** A Polyhedron Iterator is used to iterate over all polyhedron in the mesh,
	from begin() to end(), just like fr STL containers.  Its
	operator* returns a index.
	*/
	class PolyhedronIterator
	{
	public:
		/** Constructor. 
		\internal
		*/
		PolyhedronIterator() : idx_(-1) 
		{
		}
		
		PolyhedronIterator(unsigned int _idx) : idx_(_idx) 
		{
		}
		~PolyhedronIterator()
		{
		}
		/// Get  index from iterator
		unsigned int & operator*()  
		{ 
			return idx_; 
		}
		/// Get index from iterator
		unsigned int* operator->() 
		{
			return &idx_; 
		}
		/// Comparison
		bool operator==(const PolyhedronIterator & _rhs) const 
		{ 
			return idx_==_rhs.idx_;
		}
		/// Comparison
		bool operator!=(const PolyhedronIterator& _rhs) const 
		{
			return !(*this == _rhs); 
		}
		/// Pre-increment
		PolyhedronIterator & operator ++ () 
		{ 
			++idx_; 
			return *this; 
		}
		/// Pre-increment
		PolyhedronIterator & operator -- () 
		{ 
			-- idx_; 
			return *this; 
		}

		/// handle
		PolyhedronHandle handle()
		{
			return PolyhedronHandle(idx_);
		}

	private:
		unsigned int idx_;
	};


}



#endif