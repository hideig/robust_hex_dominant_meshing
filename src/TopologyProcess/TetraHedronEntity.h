#ifndef TETRA_HEDRON_ENTITY_
#define TETRA_HEDRON_ENTITY_
#include <VolumeMesh/Mesh/TetraMesh.h>
namespace TetraEntity
{

	class TetraHedronEntity
	{
	public:
		TetraHedronEntity(){}
		TetraHedronEntity(VolumeMesh::Point p_[4])
		{
			for (int i = 0; i < 4; i ++)
			{
				point[i] = p_[i];
			}
		}
		~TetraHedronEntity(){}
	public:
		//set and get function
		void setVertice(unsigned int v_[4])
		{
			for (int i = 0; i < 4; i ++)
			{
				vertice[i] = v_[i];
			}
		}
		void setPoint(VolumeMesh::Point p_[4])
		{
			for (int i = 0; i < 4; i ++)
			{
				point[i] = p_[i];
			}
		}
		void getPoint(VolumeMesh::Point *p_)
		{
			if (!p_)
			{
				p_ = new VolumeMesh::Point[4];
			}
			for (int i = 0; i < 4; i ++)
			{
				p_[i] = point[i];
			}
		}
	public:
	private:
		unsigned int vertice[4];
		VolumeMesh::Point point[4];
	};
}
#endif
/*
face is named by the opposite point
         0
        /|\
	   / | \
	  /  |  \
	 /   |   \ 
  1 /_ _ |_ _ \ 2
	\    |    /
     \   |   /
	  \  |  /
	   \ | /
	    \|/
		 3
*/   