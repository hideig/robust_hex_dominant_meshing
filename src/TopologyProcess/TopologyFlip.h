#ifndef TOPOLOGY_FLIP_
#define TOPOLOGY_FLIP_
#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TopologyProcess/TetraHedronEntity.h>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>
#include <fstream>
#include <time.h>
#include <Windows.h>
using namespace TetraEntity;
using namespace std;
class TopologyFlip
{
public:
	TopologyFlip()
	{
		isOutPut = true;
	}
	TopologyFlip(VolumeMesh::TetraMesh * tetraMesh_)
	{
		tetraMesh = tetraMesh_;
		isOutPut = true;
	}
	~TopologyFlip(){}
public:
	//set and get function
	void setTetraMesh(VolumeMesh::TetraMesh * tetraMesh_, std::string entityName_)
	{
		tetraMesh = tetraMesh_;
		entityName = entityName_;
	}

	std::set<VolumeMesh::TetraHandle> getTetrasNeedToFlip()
	{
		return tetras_NeedToFlip;
	}

	void setIsOutPut(bool isOutPut_)
	{
		isOutPut = isOutPut_;
	}
public:
	void tetraHedronFlip(int flipType, int attributeType, int qualityMeasureType);
	void flip_22_23();
	void flip_32();
	void flip_44();

//private:
	void flip_23(VolumeMesh::TetraMesh::HedronHandle flipTetra[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2]);
	double volume_length(VolumeMesh::Point v[4]);  //Volume-length measure
	double minimum_sine(VolumeMesh::Point v[4]);  // minimum sine measure
	double dihedral_sine(VolumeMesh::Point v[4], int k, int l); // calculate dihedral sine
	double minimun_jacobian(VolumeMesh::Point v[4]); // jacobian value measure of a tetrahedron
	double jacobian_value(VolumeMesh::Point v[4]);  // a single corner's jacobian value
	double edge_area(VolumeMesh::Point len_p1, VolumeMesh::Point len_p2, std::vector<VolumeMesh::Point> polygon);// ||p1 - p2|| / |area|
	double area_length(VolumeMesh::Point triangle[3]);  // |area|/(a2 + b2 + c2)
	//check if the quality of the tetra1 is better than tetra2
	bool is_optimized(std::vector<TetraEntity::TetraHedronEntity> tetra1, std::vector<TetraEntity::TetraHedronEntity> tetra2); 
	//check if the combination of the tetras has been flipped
	bool is_tetraComb_flipped(std::vector<VolumeMesh::TetraMesh::HedronHandle> tetra, std::string &tetraComb);
	// check if the combination of two tetrahedron is convex
	bool is_convex(VolumeMesh::TetraMesh::HalfFaceHandle oneOfsharedHalfFace);  

	void dataSave_flip22(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra, 
		                 VolumeMesh::TetraMesh::HedronHandle flipTetraH[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2]);
	void dataSave_flip23(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra, 
		                 VolumeMesh::TetraMesh::HedronHandle flipTetraH[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2]);
	void dataSave_flip32(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra,
		                 VolumeMesh::Point p[5], std::vector<VolumeMesh::TetraMesh::HedronHandle> edgeStar);

	void oldTetraQDS();
	void newTetraQDS();
	void updateTetraQuality(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra,
		                    std::vector<VolumeMesh::TetraMesh::HedronHandle> flipTetraHV);
	double getQualityValue(VolumeMesh::Point p[4]);
	void QualityDistrbtStatistic(double qVal, int * arrayQD);
private:
	double THRESHOLD;  //阈值
	//double duration;   //耗时
	//clock_t start, flinish;

	LARGE_INTEGER li_start, li_finish;
	double PCFreq;
	double timeCounter;

	int flipType;  //1: flip22  2:flip23  3:flip32  4:flip44
	int attributeType;  // 0x0001_volumeLength  0x0002_areaLength  0x0004_edge_area_ratio  0x0008_edge_ratio
	                    // 0x0010_corner_jac  0x0020_tetra_jac  0x0040_dihedral_sin
	int qualityMeasureType;   //  0x0001_volumeLength   0x0010_minimunSin   0x0100_jacobian

	std::string entityName;
	VolumeMesh::TetraMesh * tetraMesh;
	
	VolumeMesh::TetraMesh::HalfEdgeHandle edgeOfCoplanarFaces;   //a half edge of shared edge of coplanar faces
	std::set<std::string> flippedTetraComb;  // recording the combinations of flipped tetrahedrons 
	std::ofstream outfile, outfile2, outfile3, outfile_tc;   //write file fstream

	int oldQD[5], newQD[5];   // 质量分布统计
	int totalTetraNum;        // 总的四面体数
	int FlipTetraPairNum;     // 做翻转操作的tetra对的数目
	int FlipTimes;            // flip翻转的次数，包括回翻的flip次数
	int FlipTimes_succeeded;  // 翻转成功的flip次数
	int FlipTimes_failed;     // 翻转失败的flip次数
	std::set<VolumeMesh::TetraHandle> tetras_NeedToFlip;
	std::vector<double> tetraQuality;

	bool isOutPut;
};
#endif