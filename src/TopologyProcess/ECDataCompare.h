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
	bool is_threshold;     //�Ƿ�Ҫ��ֵ�Ƚ�

	int charaNum;                     //svm������
	std::vector<int> ECFlag;          //��������־
	std::vector<double> ECCharacter;  //svm����
	std::vector<int> predictVec;      //Ԥ���־
	std::vector<double> predictRate;  //Ԥ��ɹ�ʧ�ܵİٷֱ�,0:i*2, 1:i*2+1,iΪpredictVec�е�����

	int ecTime;           // �ܵı���������
	int succ_time;        // �������ɹ��Ĵ���
	int do_predict_not;   // Ӧ�ñ���������svmԤ��Ϊ�����Ĵ���
	int not_predict_do;   // ��Ӧ����������svmԤ��Ϊ�������Ĵ���
	int do_predict_do;
	int not_predict_not;

	std::ofstream outfile;
};
#endif