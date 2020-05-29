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

#ifndef _VOLUME_MESH_MESHIO_H_
#define _VOLUME_MESH_MESHIO_H_

//---------------------------------------------------------------------------------------------------------------------

#include <fstream>
#include <VolumeMesh/Mesh/BaseMesh.h>
#include <VolumeMesh/Mesh/TetraMesh.h>
#include <VolumeMesh/Mesh/HexadMesh.h>

/// namespace VolumeMesh
namespace VolumeMesh
{
//---------------------------------------------------------------------------------------------------------------------

/// namespace IO
namespace IO
{
//---------------------------------------------------------------------------------------------------------------------

	//bool read_mesh2Tet(TetraMesh &mesh, std::string fileName);
	std::string readline(std::ifstream * infile, int * linenumber);
	std::string find_next_sub_str(std::string & str);
	std::string trim_str(std::string & str);
	bool read_mesh(BaseMesh * & mesh, std::string fileName);
	bool read_medit_mesh(BaseMesh * & mesh, std::string fileName);
	bool read_hex(BaseMesh * & mesh, std::string fileName);
	bool read_mesh(TetraMesh & mesh, std::string fileName);
	bool read_mesh(HexadMesh & mesh, std::string fileName);
	bool write_medit_mesh(BaseMesh * mesh, const std::string & fileName);

//---------------------------------------------------------------------------------------------------------------------
}	// namespace IO
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolumeMesh
//---------------------------------------------------------------------------------------------------------------------

#endif