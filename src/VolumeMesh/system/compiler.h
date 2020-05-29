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

#ifndef _VOLUME_MESH_COMPILER_H_
#define _VOLUME_MESH_COMPILER_H_

//---------------------------------------------------------------------------------------------------------------------

#if defined(ACGMAKE_STATIC_BUILD)
#  define VM_STATIC_BUILD 1
#endif

//---------------------------------------------------------------------------------------------------------------------

#if defined(_DEBUG) || defined(DEBUG)
#  define VM_DEBUG
#endif

//---------------------------------------------------------------------------------------------------------------------

// Workaround for Intel Compiler with MS VC++ 6
#if defined(_MSC_VER) && \
	( defined(__ICL) || defined(__INTEL_COMPILER) || defined(__ICC) )
#  if !defined(__INTEL_COMPILER)
#    define __INTEL_COMPILER __ICL
#  endif
#  define VM_USE_INTEL_COMPILER 1
#endif

// ---------------------------------------------------------MS Visual C++----------------------------------------------
// Compiler _MSC_VER
// .NET 2002 1300 
// .NET 2003 1310
// .NET 2005 1400
#if defined(_MSC_VER) && !defined(VM_USE_INTEL_COMPILER)
#  if (_MSC_VER == 1300)
#    define VM_CC_MSVC
#    define VM_TYPENAME
#    define VM_OUT_OF_CLASS_TEMPLATE       0
#    define VM_PARTIAL_SPECIALIZATION      0
#    define VM_INCLUDE_TEMPLATES           1
#  elif (_MSC_VER == 1310)
#    define VM_CC_MSVC
#    define VM_TYPENAME
#    define VM_OUT_OF_CLASS_TEMPLATE       1
#    define VM_PARTIAL_SPECIALIZATION      1
#    define VM_INCLUDE_TEMPLATES           1
#    define VM_DEVELOPER_CHUHUXIAN         1
#    define VM_DEVELOPER_XIAOSHENCHEN      1
#  elif (_MSC_VER >= 1400) // settings for .NET 2005 (NOTE: not fully tested)
#    pragma warning(disable : 4996)
#    define VM_TYPENAME
#    define VM_OUT_OF_CLASS_TEMPLATE       1
#    define VM_PARTIAL_SPECIALIZATION      1
#    define VM_INCLUDE_TEMPLATES           1
#  else
#    error "Version 7 (.NET 2002) or higher of the MS VC++ is required!"
#  endif
//   currently no windows dll supported
#  define VM_STATIC_BUILD 1
#  if defined(_MT)
#    define VM_REENTRANT 1
#  endif
#  define VM_CC "MSVC++"
#  define VM_CC_VERSION _MSC_VER
// Does not work stable because the define _CPPRTTI sometimes does not exist,
// though the option /GR is set!? 
#  if defined(__cplusplus) && !defined(_CPPRTTI)
#    error "Enable Runtime Type Information (Compiler Option /GR)!"
#  endif
//#  if !defined(_USE_MATH_DEFINES_)
//#    error "You have to define _USE_MATH_DEFINES_ in the compiler settings!"
//#  endif
// ------------------------------------------------------------- Borland C --------------------------------------------
#elif defined(__BORLANDC__)
#  error "Borland Compiler are not supported yet!"
// ------------------------------------------------------------- GNU C/C++ --------------------------------------------
#elif defined(__GNUC__) && !defined(__ICC)
#  define VM_CC_GCC
#  define VM_GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 )
#  define VM_GCC_MAJOR                __GNUC__
#  define VM_GCC_MINOR                __GNUC_MINOR__
#  if (VM_GCC_VERSION >= 30200)
#    define VM_TYPENAME typename
#    define VM_OUT_OF_CLASS_TEMPLATE  1
#    define VM_PARTIAL_SPECIALIZATION 1
#    define VM_INCLUDE_TEMPLATES      1
#  else
#    error "Version 3.2.0 or better of the GNU Compiler is required!"
#  endif
#  if defined(_REENTRANT)
#    define VM_REENTRANT 1
#  endif
#  define VM_CC "GCC"
#  define VM_CC_VERSION VM_GCC_VERSION
// ------------------------------------------------------------- Intel icc --------------------------------------------
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#  define VM_CC_ICC
#  define VM_TYPENAME typename
#  define VM_OUT_OF_CLASS_TEMPLATE  1
#  define VM_PARTIAL_SPECIALIZATION 1
#  define VM_INCLUDE_TEMPLATES      1
#  if defined(_REENTRANT) || defined(_MT)
#    define VM_REENTRANT 1
#  endif
#  define VM_CC "ICC"
#  define VM_CC_VERSION __INTEL_COMPILER
//   currently no windows dll supported
#  if defined(_MSC_VER) || defined(WIN32)
#    define VM_STATIC_BUILD 1
#  endif
// ------------------------------------------------------ MIPSpro Compiler --------------------------------------------
#elif defined(__MIPS_ISA) || defined(__mips)
// _MIPS_ISA                    
// _COMPILER_VERSION            e.g. 730, 7 major, 3 minor
// _MIPS_FPSET                  32|64
// _MIPS_SZINT                  32|64
// _MIPS_SZLONG                 32|64
// _MIPS_SZPTR                  32|64
#  define VM_CC_MIPS
#  define VM_TYPENAME typename
#  define VM_OUT_OF_CLASS_TEMPLATE    1
#  define VM_PARTIAL_SPECIALIZATION   1
#  define VM_INCLUDE_TEMPLATES        0
#  define VM_CC "MIPS"
#  define VM_CC_VERSION _COMPILER_VERSION
// ------------------------------------------------------------------ other unknow compiler ---------------------------
#else
#  error "You're using an unsupported compiler!"
#endif

//---------------------------------------------------------------------------------------------------------------------
#endif // _VOLUME_MESH_COMPILER_H_ defined

