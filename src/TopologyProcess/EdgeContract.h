#ifndef TOPOLOGY_PROCESS_EDGE_CONTRACT
#define TOPOLOGY_PROCESS_EDGE_CONTRACT

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TopologyProcess/TetraHedronEntity.h>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>
#include <fstream>
#include <time.h>
#include <Windows.h>

namespace EdgeContractProcess
{
	enum {Q_VOL_LEN=0, Q_JAC_VAL, Q_MIN_SIN, Q_SQUARE_ROOT};
	class EdgeTH
	{
	public:
		EdgeTH(VolumeMesh::TetraMesh::HalfEdgeHandle he_)
		{
			heh_ = he_;
		}
		~EdgeTH(){}

		VolumeMesh::TetraMesh::HalfEdgeHandle half_edge_handle()
		{
			return heh_;
		}

		bool is_succeed()
		{
			return isSuc;
		}
		void set_suc_flag(bool isSuc_)
		{
			isSuc = isSuc_;
		}
	protected:
	private:
		VolumeMesh::TetraMesh::HalfEdgeHandle heh_;
		bool isSuc;
	};

	class EdgeContract
	{
	public:
		EdgeContract()
		{
			mesh = NULL;
		}
		EdgeContract(VolumeMesh::TetraMesh * mesh_, std::string entityName_)
		{
			mesh = mesh_;
			//undo_mesh = *mesh_;
			entityName = entityName_;
		}
		~EdgeContract(){}

		void setMesh(VolumeMesh::TetraMesh * mesh_, std::string entityName_)
		{
			mesh = mesh_;
			//undo_mesh = *mesh_;
			entityName = entityName_;
		}

		VolumeMesh::TetraMesh getMesh()
		{
			return * mesh;
		}

		void undo()
		{
			//*mesh = undo_mesh;
		}

		std::vector<EdgeTH> getContractEdge()
		{
			return edgeVec;
		}

		void initialize();
		void edgeContractData();

		//-------------------hedron quality calculate-----------------------//
		double volume_length(VolumeMesh::Point v[4]);  //Volume-length measure
		double minimum_sine(VolumeMesh::Point v[4]);  // minimum sine measure
		double minimum_sine_new(VolumeMesh::Point v[4]);  // minimum sine measure
		double square_root(VolumeMesh::Point v[4]);  // square root of radius ratio
		double dihedral_sine(VolumeMesh::Point v[4], int k, int l); // calculate dihedral sine
		double dihedral_sine(VolumeMesh::TetraMesh::HalfEdgeHandle heh_); // calculate dihedral sine
		double minimun_jacobian(VolumeMesh::Point v[4]); // jacobian value measure of a tetrahedron
		double jacobian_value(VolumeMesh::Point v[4], int vIdx);  // a single corner's jacobian value
		void jac_val_around_edge(VolumeMesh::TetraMesh::HalfEdgeHandle heh, std::vector<double> &jacValVec);  //收缩边一周的四面体网格的雅克比值计算
		double hedron_quality(VolumeMesh::Point p[4], int Q_Type = Q_VOL_LEN);
		void best_worst_quality(std::vector<VolumeMesh::TetraMesh::HedronHandle> tetraVec, 
			                    double &qVal_best, double &qVal_worst, int Q_Type = Q_VOL_LEN);
		// ratio of shortest side and longest side
		double side_ratio(VolumeMesh::TetraHandle handle);

	protected:
	private:
		std::string entityName;
		VolumeMesh::TetraMesh * mesh;
		//VolumeMesh::TetraMesh undo_mesh;
		std::vector<EdgeTH> edgeVec;
		std::set<std::string> _ephc;

		std::ofstream outfile, outfile1, outfile2, outfile_time;

		//double duration, duration_undo, duration_extra, duration_quality;   //耗时
		//clock_t start, finish, start_t, finish_t, start_q, finish_q;

		// 计时变量
		LARGE_INTEGER li_start, li_finish, li_start_t, li_finish_t, li_start_q, li_finish_q;
		double PCFreq;
		double c_ec, c_ec_undo, c_ec_extra, c_ec_quality;
	};
}

#endif