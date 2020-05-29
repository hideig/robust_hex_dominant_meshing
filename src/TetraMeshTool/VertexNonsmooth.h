#ifndef TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_
#define TETRA_MESH_TOOL_VERTEX_NONSMOOTHING_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/NonsmoothTypes.h>
#include <TetraMeshTool/NonsmoothQuality.h>
#include <TetraMeshTool/NonsmoothMath.h>
#include <map>
#include <set>
namespace VolumeMesh
{
	//define input points
	//define the points' topological structure
	class Nonsmoother
	{
	public:
		Nonsmoother(){}
		Nonsmoother(TetraMesh *_tmesh):tmesh_(_tmesh){}
		Nonsmoother(NonsmoothPoint *_points, IncidentPoints *_incident_points, int _point_num):points_(_points), vertex_num_(_point_num){}
		Nonsmoother(TetraMesh * _tmesh, int _quality_kind, float _thredhold, int _smoothing_type):tmesh_(_tmesh), quality_kinds_(_quality_kind), 
			thredhold_(_thredhold), smoothing_type_(_smoothing_type){}

		~Nonsmoother(){
			//delete points_;
			//delete old_points_;
			////delete incident_points_;
			//for(int i = 0; i < vertex_num_; ++ i)
			//{
			//	delete [] incident_tetradrons2_[i].incident_tetrahedrons;
			//}
			//delete incident_tetradrons2_;
			//delete tetrahedrons2_;
		}
		void freePointer()
		{
			delete points_;
			delete old_points_;
			//delete incident_points_;
			for(int i = 0; i < vertex_num_; ++ i)
			{
				delete [] incident_tetradrons2_[i].incident_tetrahedrons;
			}
			delete incident_tetradrons2_;
			delete tetrahedrons2_;
		}
		void setTmesh(TetraMesh * _tmesh)
		{
			tmesh_ = _tmesh;
		}
		void setPoints(NonsmoothPoint *_points)
		{
			points_ = _points;
		}
		void upDate()
		{
			TetraMesh::PointIter p_iter;
			TetraMesh::PointHandle ph;
			int index;
			TetraMesh::Point point;
			for (p_iter = tmesh_->points_begin(); p_iter != tmesh_->points_end(); ++ p_iter)
			{
				ph = p_iter.handle();
				index = ph.idx();
				point[0] = points_[index].point_coord[0];
				point[1] = points_[index].point_coord[1];
				point[2] = points_[index].point_coord[2];
				tmesh_->set_point(ph, point);
			}
		}
		void setPointNum(int _point_num)
		{
			vertex_num_ = _point_num;
		}
		void setQualityKind(int _quality_kind)
		{
			quality_kinds_ = _quality_kind;
		}
		void setSmoothingType(int _smoothing_type)
		{
			smoothing_type_ = _smoothing_type;
		}
		void setThreadhold(float _threadhold)
		{
			thredhold_ = _threadhold;
		}
		void setBoundaryChange(bool  _boundary)
		{
			boundaryChange_ = _boundary;
		}
		NonsmoothPoint * getNewPoints()
		{
			return points_;
		}
		NonsmoothPoint * getOldPoints()
		{
			return old_points_;
		}
		int getMaxIncidentNum()
		{
			//记录每个点相邻四面体最多的个数
			max_incident_num = 0;
			for(int i = 0; i < vertex_num_; ++ i)
			{
				if(incident_tetradrons2_[i].incident_hedron_num > max_incident_num)
					max_incident_num = incident_tetradrons2_[i].incident_hedron_num;
			}	
			return max_incident_num;
		}
		int getMaxIncidentVertex()
		{
			int incident_points[MAXINCIDENTPOINT];
			int incident_hedron_num;
			max_incident_vertex = 0;
			int incident_point_num;
			int point_index;
			int hedron_index;
			for (int l = 0; l < vertex_num_; ++ l)
			{
				incident_hedron_num = incident_tetradrons2_[l].incident_hedron_num;
				incident_point_num = 0;
				for (int i = 0; i < incident_hedron_num; ++ i)
				{
					hedron_index = incident_tetradrons2_[l].incident_tetrahedrons[i];
					for (int j = 0; j < 4; ++ j)
					{
						point_index = tetrahedrons2_[hedron_index].four_points_index[j];
						if (point_index != l)
						{
							int k;
							for (k = 0; k < incident_point_num; ++ k)
							{
								if (incident_points[k] == point_index)
									break;
							}
							if(k > incident_point_num - 1)
							{
								incident_points[incident_point_num] = point_index;
								++ incident_point_num;
							}
						}
					}
				}
				if (incident_point_num > max_incident_vertex)
					max_incident_vertex = incident_point_num;
			}
			return max_incident_vertex;
		}
		float getOldQuality()
		{
			return old_quality_;
		}
		float getNewQuality()
		{
			return new_quality_;
		}
		int getSmoothTime()
		{
			return smooth_time_;
		}
		int getPointNum()
		{
			return vertex_num_;
		}
		int getHedronNum()
		{
			return hedron_num_;
		}
		IncidentTetrahedrons2 * getIncidentHedron()
		{
			return incident_tetradrons2_;
		}
		OneTetrahedron2 * getHedron()
		{
			return tetrahedrons2_;
		}
		int * getPointcolor()
		{
			return point_color_;
		}
		int getGroupNum()
		{
			return group_num_;
		}
		void deleteSpace()
		{
			delete points_;
			delete old_points_;
			delete tetrahedrons2_;
			delete incident_tetradrons2_;
		}
		void setValue(int _vertex_num, int _hedron_num, double * _points, int * _hedron, int * _incident_hedron, int * _incident_hedron_num, int * _incident_point, int * _incident_point_num);
		
		void getActiveSet(int _vertex_index, double _worst_quality, NonsmoothPoint *_active_gradient, int & _active_num, HedronGradient * _gradients);
		void findDirection(int _vertex_index, NonsmoothPoint * _active_gradient, int _active_num, NonsmoothPoint &_direciton);
		void findBasis(NonsmoothPoint * _M, int _M_num, NonsmoothPoint * _S, int _S_num, NonsmoothPoint * _B, int &_B_num);
		void RecylefindBasis(NonsmoothPoint * _M, int _M_num, NonsmoothPoint * _S, int _S_num, NonsmoothPoint * _B, int &_B_num);
		double getInitialAlpha(int _vertex_index, NonsmoothPoint _direciton, double _r, double _worst_quality, HedronGradient * _gradients);
		void saveData(std::string  _file_name);
		void readData(std::string _file_name);
		double getMinQuality();
		double getMinQuality(int _vertex_index);
		void initialinformation();
		void intiialoneinformation(int _vertex_num);
		void getGradient(int  vertex_index, int hedron_index, HedronGradient &hedron_gradient_);
		void nonsmoothLineSearch(int _vertex_index, NonsmoothPoint _direction, double _r, double &_alpha, double _worst_quality);
		bool NonsmoothVertice(int _vertex_index, double &_out_worst_quality);
		bool NonsmoothVertice(double p[3], double incidenttet[][9], int tetctn, double newp[3]);
		bool NonsmoothVertex (double thredhold);
		NonsmoothPoint getFacenormal(NonsmoothPoint p1, NonsmoothPoint p2, NonsmoothPoint p3, NonsmoothPoint p4);
		bool laplacianSmooth(int _vertex_index, double & worst_quality, double thredholds);
		void pointGroup();




	private:
		TetraMesh *tmesh_;
		NonsmoothPoint *points_;
		NonsmoothPoint *old_points_;
		OneTetrahedron2 *tetrahedrons2_;

		IncidentTetrahedrons2 *incident_tetradrons2_;
		NonsmoothQuality nonsmoothquality_;
		bool boundaryChange_;

		int * point_color_;
		int quality_kinds_;
		float old_quality_;
		float new_quality_;
		float thredhold_;
		int smoothing_type_;
		
		int vertex_num_;
		int hedron_num_;
		int max_incident_num;
		int max_incident_vertex;
		int smooth_time_;
		int group_num_;
	};
}
#endif