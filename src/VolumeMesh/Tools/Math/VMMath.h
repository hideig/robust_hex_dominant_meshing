/*---------------------------------------------------------------------------------------------------------------------
*  _               _  _             _   *
*   \             /  | \           / |  *     VolumeMeshTools : The Open source for providing some operations on 
*    \           /   |  \         /  |  *                       VolumeMesh 
*  	  \         /    |   \       /   |  *     Copyright(C) 2010 by Computer Aided Designed Group
*	   \       /     |    \     /    |  *     State Key Lab. of CAD & CG, Zhejiang University
*	    \     /      |     \   /     |  *
*	     \   /       |      \_/      |  *
*		  \_/        |               |_ * 
*                                       *
*-----------------------------------------------------------------------------------------------------------------------
* License
*
*    This file is part of VolumeMeshTools.
*
*    VolumeMeshTools is free software: you can redistribute it and/or modify       
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
*
/---------------------------------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------------------------------
*  File Information                                                                                                   *
*  written by Chuhua Xian                                                                                             *
*  Email : chuhuaxian@gmail.com                                                                                       *
*  modified day : 2010/7/28.                                                                                          *
*--------------------------------------------------------------------------------------------------------------------*/


#ifndef _VOLUME_MESH_MATH_
#define _VOLUME_MESH_MATH_


#define VMT_EPS 10e-10
#define VMT_MAX 10e20

#include <math.h>


namespace VolumeMesh
{
	namespace VMMath
	{	
//---------------------------------------------------------------------------------------------------------------------		

		const static double VM_EPS = 1e-10;
		const static double VM_MAX = 1e20;
		const static double VM_PI  = 3.14159265358979323846;

		template <typename T>
		inline T zeroValue()
		{
			return T(0.0);
		}
		/**
		* this namespace provides some operations of linear algebra algorithm, all functions are written in 
		* template form for generic applications
		*/
		namespace LinearAgebra
		{
			//---------------------------------------------------------------------------------------------------------------------
			/** calculate the determination of a 3 * 3 matrix
			*   \param in _A[3] the matrix is stored in a 3 - vector
			*   \return the determination 
			*/
			template <typename T>
			inline T determination3x3(const T _A[3][3])
			{
				return _A[0][0] * _A[1][1] * _A[2][2] + 
					   _A[1][0] * _A[2][1] * _A[0][2] +
					   _A[0][1] * _A[1][2] * _A[2][0] -
					   _A[0][2] * _A[1][1] * _A[2][0] -
					   _A[0][1] * _A[1][0] * _A[2][2] -
					   _A[1][2] * _A[2][1] * _A[0][0];

			}

			/** 
			* solve three order equation Ax = b using the crame's law
			* @param [in] : A[3][3]. the matrix of Ax = b
			* @param [in] : b[3]. the b of Ax = b
			* @return : bool. true for success, false for not
			**/
			template <typename T>
			inline bool solveLinearSystem3x3(T _A[3][3], T * _b)
			{
				T crame;
				crame = determination3x3(_A);
				if (crame == zeroValue<T>())
				{
					return false;
				}
				T A0[3][3];
				T A1[3][3];
				T A2[3][3];
				for (unsigned int i = 0; i < 3; i ++)
				{
					for (unsigned int k = 0; k < 3; k ++)
					{
						A0[i][k] = _A[i][k];
						A1[i][k] = _A[i][k];
						A2[i][k] = _A[i][k];
					}
					A0[i][0] = _b[i];
					A1[i][1] = _b[i];
					A2[i][2] = _b[i];
				}
				_b[0] = determination3x3(A0) / crame;
				_b[1] = determination3x3(A1) / crame;
				_b[2] = determination3x3(A2) / crame;

				return true;
			}

			/** 
			* linear transformation T = Mx
			**/
			template <typename Vec>
			inline Vec linearTransformation(Vec M[], Vec x)
			{
				Vec T;
				T = zeroValue<Vec>();
				unsigned int d;
				d = Vec::dim();
				for (unsigned int i = 0; i < d; i ++)
				{
					for (unsigned int k = 0; k < d; k ++)
					{
						T[i] += M[i][k] * x[k];
					}
				}
				return T;
			}

			/** 
			* linear transformation T = Mx
			**/
			template <typename T, typename Vec>
			inline void linearTransformation(const T M[3][3], Vec & x)
			{
				Vec T;
				T = zeroValue<Vec>();
				unsigned int d;
				d = Vec::dim();
				for (unsigned int i = 0; i < d; i ++)
				{
					for (unsigned int k = 0; k < d; k ++)
					{
						T[i] += M[i][k] * x[k];
					}
				}
				x = T;
			}

			

			//---------------------------------------------------------------------------------------------------------------------
		}

//---------------------------------------------------------------------------------------------------------------------

	}

}
#endif