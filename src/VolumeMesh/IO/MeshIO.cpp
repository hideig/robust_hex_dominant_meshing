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

#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <VolumeMesh/io/MeshIO.h>
#include <VolumeMesh/Mesh/TetraMesh.h>
#include <sstream>
//#include <VolumeMesh/Mesh2Tet/Mesh2Tet.h>

#define INPUTLINESIZE 2048


//using namespace std;


namespace VolumeMesh	
{
//---------------------------------------------------------------------------------------------------------------------
namespace IO	
{
	
	/** readline()   Read a nonempty line from a file.                            
	*                                                                           
	* A line is considered "nonempty" if it contains something more than white  
	* spaces.  If a line is considered empty, it will be dropped and the next   
	* line will be read, this process ends until reaching the end-of-file or a  
	* non-empty line.  Return NULL if it is the end-of-file, otherwise, return  
	* a pointer to the first non-whitespace character of the line.              
	*/

	std::string readline(std::ifstream * infile, int * linenumber)
	{
		char line[INPUTLINESIZE]; //
		//std::string line;
		char * result;
		//bool f;
		

		// Search for a non-empty line.
		do 
		{
			//result = fgets(str, INPUTLINESIZE - 1, infile);
			if (infile->getline(line, INPUTLINESIZE - 1))
			{
				if (linenumber) 
				{
					(*linenumber) ++;
				}

				result = line;
			    // Skip white spaces.
				while ((*result == ' ') || (*result == '\t')) 
				{
					result ++;
				}
			}
			else  
			{
				return std::string("");
			}

			// If it's end of line, read another line and try again.
		} while (*result == '\0');
		
		return std::string(result);
	}

	std::string trim_str(std::string & str)
	{

		int i;
		int k;
		i = 0;

		// skip the space char
		i = 0;
		k = str.size() - 1;
		while ((str[i] == ' ') || (str[i] == '\t'))
		{
			i ++;
		}

		while ((str[k] == ' ') || (str[k] == '\t'))
		{
			k --;
		}
		str = str.substr(i, k - i + 1);

		return str;
	}
	
	
	                                                                           
	/** findnextnumber()   Find the next field of a number string.              
	*                                                                          
	* Jumps past the current field by searching for whitespace or a comma, then
	* jumps past the whitespace or the comma to find the next field that looks 
	* like a number.                                                           
	*                                                                          
	*/

	std::string find_next_sub_str(std::string & str)
	{
		
		int i;
		int k;
		i = 0;

		// skip the space char
		i = 0;
		k = str.size() - 1;
		while ((str[i] == ' ') || (str[i] == '\t'))
		{
			i ++;
		}

		while ((str[k] == ' ') || (str[k] == '\t'))
		{
			k --;
		}
		str = str.substr(i, k - i + 1);

		//result = str;
		// Skip the current field.  Stop upon reaching whitespace or a comma.
		while ((str[i] != '\0') && (str[i] != '#') && (str[i] != ' ') && 
			   (str[i] != '\t') && (str[i] != ',')) 
		{
			i ++;
		}

		// Now skip the whitespace and anything else that doesn't look like a
		//   number, a comment, or the end of a line. 
		while ((str[i] != '\0') && (str[i] != '#') 
			   && (str[i] != '.') && (str[i] != '+') && (str[i] != '-')
			   && ((str[i] < '0') || (str[i] > '9'))) 
		{
				i ++;
		}

		k = i;
		while ((str[k] != '\0') && (str[k] != '#') 
			   && ((str[k] == '.') || (str[k] == '+') || (str[k] == '-')
			   ||  (str[k] == 'e') || (str[k] == 'E')
			   ||  ((str[k] >= '0') && (str[k] <= '9')))) 
		{
			k ++;
		}


		

		// Check for a comment (prefixed with `#').
		if (str[i] == '#') 
		{
			return std::string("");
		}
		std::string rlt;
		rlt = str.substr(i - 1, k - i + 1);
		str = str.substr(i - 1);
		return rlt;
	}



	//-----------------------------------------------------------------------------------------------------------------

	/** read the mesh file, currently support file format : .mesh, .tet, .hex
	*/
	bool read_mesh(BaseMesh * & mesh, std::string fileName)
	{
		if (mesh)
		{
			delete mesh;
			mesh = NULL;
		}
		
		if (fileName.find(".mesh") == fileName.size() - 5)
		{
			return read_medit_mesh(mesh, fileName);
		}
		else if (fileName.find(".hex") == fileName.size() - 4)
		{
			return read_hex(mesh, fileName);
		}
		
		return false;
	}

	/** read the .mesh file which can be open by the free software MEdit
	*/
	bool read_medit_mesh(BaseMesh * & mesh, std::string fileName)
	{
		//ensure the file format
		if ( fileName.find(".mesh") != fileName.size() - 5 )
		{
			std::cerr << "Error: File format invalid" << std::endl;
			return false;
		}

		std::ifstream fs(fileName.c_str());
		if ( !fs.is_open() )
		{
			std::cerr << "Error: File not found!" << std::endl;
			return false;
		}

		//mesh.release_mesh();

		int lineCount;
		std::string buffer;
		char * str;
		int dimension;
		unsigned int nv;
		unsigned int nf;
		unsigned int nt;
		int corners;
		int mt;
		int firstLineNumber;
		lineCount = 0;
		dimension = 0;
		nv = 0;
		nf = 0;
		nt = 0;
		firstLineNumber = 1;
		//buffer = NULL;
		str = NULL;
		PContainer pc;
		VContainer vc;
		double v[4];
		int idx[9];

		while (!(buffer = readline(&fs, &lineCount)).empty()) 
		{
			if (buffer[0] != '#') 
			{
				if (!dimension)
				{
					// Find if it is the keyword "Dimension".					
					if ((buffer.find("Dimension") == 0) || (buffer.find("DIMENSION") == 0))
					{
						// Read the dimensions

						trim_str(buffer);

						if (buffer.size() != std::string("Dimension").size())
						{
							buffer = find_next_sub_str(buffer);
							std::istringstream ss(buffer);
							ss >> dimension;
						}
						else
						{
							fs >> dimension;

						}
						if (dimension != 2 && dimension != 3) 
						{
							std::cerr << "Unknown dimension in file on line"  << lineCount << " in file \n" << fileName << std::endl;
							fs.close();
							return false;
						}
					}
				} // end of dimension

				if (dimension && (!nv))
				{
					// find the key work "Vertices"
					if (buffer == "Vertices")
					{
						fs >> nv;
						if (!nv)
						{
							return false;
						}	

						pc.resize(nv, Point(0, 0, 0));
						for (unsigned int i = 0; i < nv; i ++)
						{
							fs >> v[0] >> v[1] >> v[2] >> v[3];
							pc[i] = Point(v[0], v[1], v[2]);
						}										
					}
				}
				if (nv && (!nf))
				{
					// find the key word "Triangles" or "Quadrilaterals"
					if ((buffer == "Triangles") || (buffer == "TRIANGLES"))
					{
						fs >> nv;
						corners = 3;
					}
					else if ((buffer == "Quadrilaterals") || (buffer == "QUADRILATERALS"))
					{
						fs >> nv;
						corners = 4;						
					}	
				}
				if (nv && (!nt))
				{
					if ((buffer == "Tetrahedra") || (buffer == "TETRAHEDRA") || 
						(buffer == "Hexahedra")  || (buffer == "HEXAHEDRA"))
					{
						fs >> nt;
						if ((buffer == "Tetrahedra") || (buffer == "TETRAHEDRA"))
						{
							mt = 4;
						}
						else if ((buffer == "Hexahedra")  || (buffer == "HEXAHEDRA"))
						{
							mt = 8;
						}
						
						
						
						if (!nt)
						{
							return false;
						}

						for (unsigned int i = 0; i < nt; i ++)
						{
							if (mt == 4)
							{
								fs >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[5];
							}
							else if (mt == 8)
							{
								fs >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >> idx[6] >> idx[7] >> idx[8];
							}

							for (int k = 0; k < mt; k ++)
							{
								vc.push_back(Vertex(i * mt + k, idx[k] - 1));
							}
						}
					}					
				}


			} // end of if '#'
		} // end of while
		fs.close();

		if (!mesh)
		{
			if (mt & 0x0004)
			{
				// tetrahedral mesh
				mesh = new TetraMesh; 
			}
			else if (mt & 0x0008)
			{
				// hexahedral mesh
				mesh = new HexadMesh;
			}
			else
			{
				return false;
			}
		}
		mesh->set_point_container(pc);
		mesh->set_vertex_container(vc);

		return mesh->build_topology();
	}

	bool read_hex(BaseMesh * & mesh, std::string fileName)
	{
		trim_str(fileName);
		if (fileName.find(".hex") != fileName.size() - 4)
		{
			std::cerr << "Error : not a .hex file!" << std::endl;
			return false;
		}
		std::ifstream fs(fileName.c_str());

		int lineCount;
		std::string buffer;
		std::string next_str;
		unsigned int nv;
		unsigned int nf;
		unsigned int nt;
		int firstLineNumber;
		lineCount = 0;
		nv = 0;
		nf = 0;
		nt = 0;
		firstLineNumber = 1;

		PContainer pc;
		VContainer vc;
		double v[3];
		int idx[9];
		
		

		while (!(buffer = readline(&fs, &lineCount)).empty()) 
		{
			//trim_str(buffer);
			if (buffer[0] != '#') 
			{	
				if ((buffer[0] == 'v') || (buffer[0] == 'V'))
				{
					for (unsigned int i = 0; i < 3; i ++)
					{
						next_str = find_next_sub_str(buffer);
						std::istringstream ss(next_str);
						ss >> v[i];
					}
					pc.push_back(Point(v[0], v[1], v[2]));
					++ nv;


				}
				
				if ((buffer[0] == 'h') || (buffer[0] == 'H'))
				{
					for (unsigned int i = 0; i < 8; i ++)
					{
						next_str = find_next_sub_str(buffer);
						std::istringstream ss(next_str);
						ss >> idx[i];
						vc.push_back(Vertex(nt * 8 + i, idx[i] - 1));
					}
					++ nt;
					//break;
				}

			} // end of if '#'
		} // end of while
		fs.close();

		if (!mesh)
		{
			mesh = new HexadMesh;
		}
		mesh->set_point_container(pc);
		mesh->set_vertex_container(vc);


		return mesh->build_topology();
	}

//---------------------------------------------------------------------------------------------------------------------


	/** mesh read function for TetraMesh
	*/
	bool read_mesh(TetraMesh &mesh, std::string fileName)
	{
		//ensure the file format
		if ( fileName.find(".mesh") != fileName.size() - 5 )
		{
			std::cerr << "Error: File format invalid" << std::endl;
			return false;
		}

		std::ifstream fs(fileName.c_str());
		if ( !fs.is_open() )
		{
			std::cerr << "Error: File not found!" << std::endl;
			return false;
		}

		mesh.release_mesh();

		int lineCount;
		std::string buffer;
		char * str;
		int dimension;
		unsigned int nv;
		unsigned int nf;
		unsigned int nt;
		int corners;
		int mt;
		int firstLineNumber;
		lineCount = 0;
		dimension = 0;
		nv = 0;
		nf = 0;
		nt = 0;
		firstLineNumber = 1;
		//buffer = NULL;
		str = NULL;
		PContainer pc;
		VContainer vc;
		double v[4];
		int idx[9];

		while (!(buffer = readline(&fs, &lineCount)).empty()) 
		{
			if (buffer[0] != '#') 
			{
				if (!dimension)
				{
					// Find if it is the keyword "Dimension".					
					if ((buffer.find("Dimension") == 0) || (buffer.find("DIMENSION") == 0))
					{
						// Read the dimensions

						trim_str(buffer);

						if (buffer.size() != std::string("Dimension").size())
						{
							buffer = find_next_sub_str(buffer);
							std::istringstream ss(buffer);
							ss >> dimension;
						}
						else
						{
							fs >> dimension;

						}
						if (dimension != 2 && dimension != 3) 
						{
							std::cerr << "Unknown dimension in file on line"  << lineCount << " in file \n" << fileName << std::endl;
							fs.close();
							return false;
						}
					}
				} // end of dimension

				if (dimension && (!nv))
				{
					// find the key work "Vertices"
					if (buffer == "Vertices")
					{
						fs >> nv;
						if (!nv)
						{
							return false;
						}	

						pc.resize(nv, Point(0, 0, 0));
						for (unsigned int i = 0; i < nv; i ++)
						{
							fs >> v[0] >> v[1] >> v[2] >> v[3];
							pc[i] = Point(v[0], v[1], v[2]);
						}										
					}
				}
				if (nv && (!nf))
				{
					// find the key word "Triangles" or "Quadrilaterals"
					if ((buffer == "Triangles") || (buffer == "TRIANGLES"))
					{
						fs >> nv;
						corners = 3;
					}
					else if ((buffer == "Quadrilaterals") || (buffer == "QUADRILATERALS"))
					{
						fs >> nv;
						corners = 4;						
					}	
				}
				if (nv && (!nt))
				{
					if ((buffer == "Tetrahedra") || (buffer == "TETRAHEDRA"))
					{
						fs >> nt;
						mt = 4;
						if (!nt)
						{
							return false;
						}

						for (unsigned int i = 0; i < nt; i ++)
						{
							if (mt == 4)
							{
								fs >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4];
							}
							else if (mt == 8)
							{
								fs >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >> idx[6] >> idx[7] >> idx[8];
							}
							
							for (int k = 0; k < mt; k ++)
							{
								vc.push_back(Vertex(i * mt + k, idx[k] - 1));
							}
						}
					}					
				}
				
				
			} // end of if '#'
		} // end of while
		fs.close();

		mesh.set_point_container(pc);
		mesh.set_vertex_container(vc);

		return mesh.build_topology();
	
	}

	/** 
	 * write the medit file
	 * \param in mesh the mesh
	 * \param in fileName the file name
	 * \return bool true for success, false for failure
	 */
	bool write_medit_mesh(BaseMesh * mesh, const std::string & fileName)
	{
		std::ofstream f(fileName.c_str());
		f << "MeshVersionFormatted 1\n";
		f << "Dimension \n 3\n";
		f << "Vertices\n";
		f << mesh->n_point() << '\n';
		
		
		Point p;
		if (mesh->mesh_type() & 0x0001)
		{
		   /** 
		   * write the tetrahedral mesh
		   */
			TetraMesh * t_mesh;
			t_mesh = (TetraMesh *)mesh;
			for (TetraMesh::PointIter p_it = t_mesh->points_begin(); p_it != t_mesh->points_end(); ++ p_it)
			{
				p = t_mesh->point(p_it);
				
				f << p[0] << "   " << p[1] << "   " << p[2] << "   " << 0 << '\n'; 
			}

			f << "Triangles\n" << 1 << '\n';

			f << "Tetrahedra\n";
			f << mesh->n_polyhedron() << '\n';



			TetraMesh::HedronIter h_it;
			TetraMesh::HedronVertexIter hv_it;
			for (h_it = t_mesh->hedrons_begin(); h_it != t_mesh->hedrons_end(); ++ h_it)
			{
				for (hv_it = t_mesh->hedron_vertex_iter(h_it); hv_it; ++ hv_it)
				{
					f << t_mesh->handle_to_entity(hv_it.handle()).point_handle().idx() + 1 << "  ";
				}
				f << 0 << '\n';				 
			}
		}
		else if (mesh->mesh_type() & 0x0002)
		{
			HexadMesh * h_mesh;
			h_mesh = (HexadMesh *)mesh;
			for (HexadMesh::PointIter p_it = h_mesh->points_begin(); p_it != h_mesh->points_end(); ++ p_it)
			{
				p = h_mesh->point(p_it);

				f << p[0] << "   " << p[1] << "   " << p[2] << "   " << 0 << '\n'; 
			}

			f << "Quadrilaterals\n" << 1 << '\n';

			f << "Hexahedra\n";
			f << mesh->n_polyhedron() << '\n';



			HexadMesh::HedronIter h_it;
			HexadMesh::HedronVertexIter hv_it;
			unsigned int hidx[8];
			for (h_it = h_mesh->hedrons_begin(); h_it != h_mesh->hedrons_end(); ++ h_it)
			{
				int i = 0;
				for (hv_it = h_mesh->hedron_vertex_iter(h_it); hv_it; ++ hv_it)
				{
					hidx[i] = h_mesh->handle_to_entity(hv_it.handle()).point_handle().idx() + 1;
					++ i;
				}
				f << hidx[0] << "   " << hidx[1] << "   " << hidx[2] << "   " << hidx[3] << "   "; 
				f << hidx[4] << "   " << hidx[5] << "   " << hidx[6] << "   " << hidx[7] << "   "; 
				f << 0 << '\n';				 
			}
		}

		f << "End";

		
		f.close();
		
		return true;
	}



//---------------------------------------------------------------------------------------------------------------------
}	// namespace IO
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
}	// namespace VolueMesh
//---------------------------------------------------------------------------------------------------------------------
