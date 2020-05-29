#ifndef TOPOLOGY_PROCESS_DATA_COMPARE
#define TOPOLOGY_PROCESS_DATA_COMPARE
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <vector>

class DataCompare
{
public:
	DataCompare(){}
	DataCompare(std::string filePath_, std::string fileName_, int flipType_, int qualityMeasureType_)
	{
		filePath = filePath_;
		fileName = fileName_;
		flipType = flipType_;
		qualityMeasureType = qualityMeasureType_;
	}
	~DataCompare(){}

	void setParameters(std::string filePath_, std::string fileName_, int flipType_, int qualityMeasureType_)
	{
		filePath = filePath_;
		fileName = fileName_;
		flipType = flipType_;
		qualityMeasureType = qualityMeasureType_;
	}
public:
	void initialize();
	void process_threshold(); // 阈值相关的数据处理
	void read_predictData(std::string fileName);  // 读取.t.predict文件
	void read_qualityData(std::string fileName);  // 读取quality文件

	void data_with_threshold(double threshold);    // 加上阈值后的数据
	void data_with_different_threshold();    // 加上不同阈值后的数据
	void QualityDistrbtStatistic(std::vector<double> qualityArray, int * arrayQD);  // 质量分布统计
	void updateTetraQuality(int idx, std::vector<double> &qualityV, int *v_flag);

	double accuracyRateCalculate(std::vector<int> standardVec, std::vector<int> testVec); // 正确率计算
	int rate_do_predict_not(std::vector<int> standardVec, std::vector<int> testVec);   //应翻转，但svm预测为不翻转的情况统计
	int rate_not_predict_do(std::vector<int> standardVec, std::vector<int> testVec);   //不应翻转，但svm预测为翻转的情况统计
protected:
private:
	std::string filePath, fileName, entityName;
	int flipType;            // flip22: 0x0001  flip23: 0x0002  flip32: 0x0004  flip44: 0x0008
	int qualityMeasureType;  // VolLen: 0x0001  MinSin: 0x0002  JacVal: 0x0004

	std::string FN_predict, FN_quality;
	std::ofstream outfile;

	std::vector<int> label_quality;  // labels in flipData
	std::vector<int> label_predict;  // labels in predictData
	std::vector<int> flipTetraIdx;   // 翻转的网格索引
	std::vector<double> quality, quality_flip, quality_svm;
	std::vector<double> quality_new;
	std::vector<double> predict_1;
	std::vector<double> predict_0;

	int radix_old, radix_new;       // radix_old  参与翻转的网格数目 radix_new  翻转生成的网格数目
	int totalFlipTime;   // 总的flip翻转次数
	int do_predict_not;  // 应翻转，但svm预测为不翻转的flip次数
	int not_predict_do;  // 不应翻转，但svm预测为翻转的flip次数
	int oldQD[5], flipQD[5], svmQD[5];      // 质量分布数组 array[0]:0~0.1  array[1]:0.1~0.3  array[2]:0.3~0.5
	                     //               array[3]:0.5~0.7  array[4]:0.7~1.0
	int *quality_flip_flag, *quality_svm_flag;
	double accuracyRate; // 正确率
};
#endif