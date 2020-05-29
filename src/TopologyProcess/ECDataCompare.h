#ifndef TOPOLOGY_PROCESS_EC_DATA_COMPARE
#define TOPOLOGY_PROCESS_EC_DATA_COMPARE
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <vector>
class ECDataCompare
{
public:
	ECDataCompare(std::string ecfile_, std::string pfile_, bool threshold_)
	{
		ecfile = ecfile_;
		pfile = pfile_;
		is_threshold = threshold_;
	}
	~ECDataCompare(){}

public:
	int initialize();
	void read_predict_file();
	void read_svm_file();
	void ec_data_compare();
	void data_analyse();
	void data_analyse_threshold(double threshold);

protected:
private:
	std::string ecfile, pfile;
	std::string path;
	std::string modelName;
	bool is_threshold;     //是否要阈值比较

	int charaNum;                     //svm特征数
	std::vector<int> ECFlag;          //边收缩标志
	std::vector<double> ECCharacter;  //svm特征
	std::vector<int> predictVec;      //预测标志
	std::vector<double> predictRate;  //预测成功失败的百分比,0:i*2, 1:i*2+1,i为predictVec中的索引

	int ecTime;           // 总的边收缩次数
	int succ_time;        // 边收缩成功的次数
	int do_predict_not;   // 应该边收缩，但svm预测为不做的次数
	int not_predict_do;   // 不应边收缩，但svm预测为边收缩的次数
	int do_predict_do;
	int not_predict_not;

	std::ofstream outfile;
};
#endif