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
	void process_threshold(); // ��ֵ��ص����ݴ���
	void read_predictData(std::string fileName);  // ��ȡ.t.predict�ļ�
	void read_qualityData(std::string fileName);  // ��ȡquality�ļ�

	void data_with_threshold(double threshold);    // ������ֵ�������
	void data_with_different_threshold();    // ���ϲ�ͬ��ֵ�������
	void QualityDistrbtStatistic(std::vector<double> qualityArray, int * arrayQD);  // �����ֲ�ͳ��
	void updateTetraQuality(int idx, std::vector<double> &qualityV, int *v_flag);

	double accuracyRateCalculate(std::vector<int> standardVec, std::vector<int> testVec); // ��ȷ�ʼ���
	int rate_do_predict_not(std::vector<int> standardVec, std::vector<int> testVec);   //Ӧ��ת����svmԤ��Ϊ����ת�����ͳ��
	int rate_not_predict_do(std::vector<int> standardVec, std::vector<int> testVec);   //��Ӧ��ת����svmԤ��Ϊ��ת�����ͳ��
protected:
private:
	std::string filePath, fileName, entityName;
	int flipType;            // flip22: 0x0001  flip23: 0x0002  flip32: 0x0004  flip44: 0x0008
	int qualityMeasureType;  // VolLen: 0x0001  MinSin: 0x0002  JacVal: 0x0004

	std::string FN_predict, FN_quality;
	std::ofstream outfile;

	std::vector<int> label_quality;  // labels in flipData
	std::vector<int> label_predict;  // labels in predictData
	std::vector<int> flipTetraIdx;   // ��ת����������
	std::vector<double> quality, quality_flip, quality_svm;
	std::vector<double> quality_new;
	std::vector<double> predict_1;
	std::vector<double> predict_0;

	int radix_old, radix_new;       // radix_old  ���뷭ת��������Ŀ radix_new  ��ת���ɵ�������Ŀ
	int totalFlipTime;   // �ܵ�flip��ת����
	int do_predict_not;  // Ӧ��ת����svmԤ��Ϊ����ת��flip����
	int not_predict_do;  // ��Ӧ��ת����svmԤ��Ϊ��ת��flip����
	int oldQD[5], flipQD[5], svmQD[5];      // �����ֲ����� array[0]:0~0.1  array[1]:0.1~0.3  array[2]:0.3~0.5
	                     //               array[3]:0.5~0.7  array[4]:0.7~1.0
	int *quality_flip_flag, *quality_svm_flag;
	double accuracyRate; // ��ȷ��
};
#endif