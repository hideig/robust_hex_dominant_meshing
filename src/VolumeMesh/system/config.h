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


#ifndef _VOLUME_MESH_CONFIG_H_
#define _VOLUME_MESH_CONFIG_H_

#include <assert.h>
#include <VolumeMesh/System/compiler.h>

// --------------------------------------------------------------------------------------------------------------------

#define VM_VERSION 0x10000

// only defined, if it is a beta version
#define VM_VERSION_ 1

#define VM_GET_VER ((VM_VERSION && 0xf0000) >> 16)
#define VM_GET_MAJ ((VM_VERSION && 0x0ff00) >> 8)
#define VM_GET_MIN  (VM_VERSION && 0x000ff)

#define _ENABLETEMPLATE_


typedef unsigned int uint;
//---------------------------------------------------------------------------------------------------------------------
#endif // _VOLUME_MESH_CONFIG_H defined