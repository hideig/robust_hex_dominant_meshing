#ifndef TOPOLOGY_PROCESS_EC_DATA_PROCESS
#define TOPOLOGY_PROCESS_EC_DATA_PROCESS
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <vector>
class ECDataProcess
{
public:
	//ECDataProcess()
	//{

	//}
	ECDataProcess(std::string fileName_, int svm_attr_, double threshold_, int svm_tetra_num_)
	{
		fileName = fileName_;
		svm_attr = svm_attr_;
		threshold = threshold_;
		svm_tetra_num = svm_tetra_num_;
	}
	~ECDataProcess()
	{

	}
	void initialize();
	void edge_contract_data_process();
	void save_svm_data();
	void save_svm_data_single_line(int succ, std::vector<double> attr, std::ofstream &outfile);
	void tetra_range_QIWT();
	void tetra_range_QIET();
	void array_idx_TetraRange_QIWT(double row, double column, int &i, int &j);
	void array_idx_TetraRange_QIET(double row, int column, int &i, int &j);
protected:
private:
	int svm_attr;
	int line_counter;
	int succ_num;
	double threshold;
	int svm_tetra_num;
	std::string fileName;
	std::vector<std::vector<double>> jacValAroudEdge;
	std::vector<int> ECFlag;
	std::vector<int> ETVec;  // pf_tn, pt_tn, ec_tn
	std::vector<double> QVec;  //qVal_best, qVal_worst, qVal_ec_best, qVal_ec_worst
	std::vector<int> locationVec;
	std::vector<double> QImproveVec;
	std::vector<double> minSineVec;
	std::vector<double> squareRoot;
	std::vector<double> sideRatio;  // sr_large	sr_small, sr_aver, sr_ec_large, sr_ec_small, sr_ec_aver, sr_ep_large, sr_ep_small, sr_ep_aver

	int TetraRange_QIWT[15][10];
	int TetraRange_QIET[9][18];

	std::ofstream outfile;
	std::ofstream ofSortTN[9];
};
#endif