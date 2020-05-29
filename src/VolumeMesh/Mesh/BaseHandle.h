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


#ifndef _VOLUME_MESH_BASE_HANDLE_H_
#define _VOLUME_MESH_BASE_HANDLE_H_
//---------------------------------------------------------------------------------------------------------------------


#include <VolumeMesh/system/config.h>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------
	class BaseHandle
	{
	public:

		explicit BaseHandle(int _idx=-1) : idx_(_idx) 
		{
		}

		// Get the underlying index of this handle
		int idx() const 
		{
			return idx_;
		}

		// The handle is valid iff the index is not equal to -1.
		bool is_valid() const 
		{
			return idx_ != -1; 
		}

		// reset handle to be invalid
		void reset() 
		{
			idx_= -1; 
		}

		// reset handle to be invalid
		void invalidate() 
		{
			idx_ = -1; 
		}

		bool operator == (const BaseHandle& _rhs) const 
		{ 
			return (idx_ == _rhs.idx_); 
		}

		bool operator != (const BaseHandle& _rhs) const 
		{ 
			return (idx_ != _rhs.idx_); 
		}

		bool operator < (const BaseHandle& _rhs) const 
		{ 
			return (idx_ < _rhs.idx_); 
		}

		operator int() const
		{
			return idx_;
		}

		// this is to be used only by the iterators
		void __increment() 
		{ 
			++ idx_; 
		}

		void __decrement() 
		{
			-- idx_; 
		}

		void setIdx(int idx)
		{
			idx_ = idx;
		}

	private:
		int idx_; 
	};

//---------------------------------------------------------------------------------------------------------------------
} // namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------

#endif // _VOLUME_MESH_BASE_HANDLE_H_ defined
