#include <TopologyProcess/DataCompare.h>
#include <math.h>
#include <iomanip>
#include <windows.h>
#include <sstream>
using namespace std;

void DataCompare::initialize()
{
	// get predict file name and quality file name
	FN_quality = filePath;
	FN_predict = filePath;
	FN_quality += "\\Q_flip";
	FN_predict += "\\svm_flip";
	entityName = "flip";

	if (flipType & 0x0001)
	{
		FN_quality += "22_";
		FN_predict += "22_";
		entityName += "22_";
	}
	else if (flipType & 0x0002)
	{
		FN_quality += "23_";
		FN_predict += "23_";
		entityName += "23_";
	}
	else if (flipType & 0x0004)
	{
		FN_quality += "32_";
		FN_predict += "32_";
		entityName += "32_";
	}
	else if (flipType & 0x0008)
	{
		FN_quality += "44_";
		FN_predict += "44_";
		entityName += "44_";
	}

	FN_quality += fileName;
	FN_predict += fileName;
	entityName += fileName;

	if (qualityMeasureType & 0x0001)
	{
		FN_quality += "_VolLen";
		FN_predict += "_VolLen";
		entityName += "_VolLen";
	}
	if (qualityMeasureType & 0x0002)
	{
		FN_quality += "_MinSin";
		FN_predict += "_MinSin";
		entityName += "_MinSin";
	}
	if (qualityMeasureType & 0x0004)
	{
		FN_quality += "_JacVal";
		FN_predict += "_JacVal";
		entityName += "_JacVal";
	}
	FN_predict += ".t.predict";
}

void DataCompare::process_threshold()
{
	initialize();
	read_predictData(FN_predict);
	read_qualityData(FN_quality);

	// check data
	if (label_predict.size() != label_quality.size())
	{
		MessageBox(0,TEXT("Wrong Data!"),TEXT("DataError"),MB_OK);
		return;
	}

	data_with_different_threshold();
}

void DataCompare::read_predictData(std::string fileName)
{
	std::ifstream in_predict;
	in_predict.open(fileName.c_str());
	if (!in_predict)
	{
		return;
	}

	std::string text;
	int label;
	double p1, p0;
	getline(in_predict, text);
	getline(in_predict, text);
	while(text != "\n" && text.size() != 0)
	{
		std::istringstream sin(text);
		sin >> label >> p0 >> p1;
		label_predict.push_back(label);
		predict_0.push_back(p0);
		predict_1.push_back(p1);
		text.clear();
		getline(in_predict, text);
	}
}

void DataCompare::read_qualityData(std::string fileName)
{
	std::ifstream in_quality;
	in_quality.open(fileName.c_str());
	if (!in_quality)
	{
		return;
	}

	std::string text, s_flipType, s_qualityType;
	int label, index, count;
	double Q_old;
	double *Q_new;
	getline(in_quality, text);
	std::istringstream sin(text);
	sin >> s_flipType >> s_qualityType;

	// get original tetra quality values
	getline(in_quality, text);
	sin.clear();
	sin.str(text);
	sin >> count;
	quality.resize(count);
	for (int i = 0; i < count; i++)
	{
		getline(in_quality, text);
		sin.clear();
		sin.str(text);
		sin >> Q_old;
		quality[i] = Q_old;
	}

	// get flip tetra index and new tetra quality values
	if (s_flipType == "FLIP22")
	{
		radix_old = radix_new = 2;
		Q_new = new double[2];
		getline(in_quality, text);
		while (text != "\n" && text.size() != 0)
		{
			sin.clear();
			sin.str(text);
			// get label
			sin >> label;
			label_quality.push_back(label);
			// get flip tetra index
			sin >> index;
			flipTetraIdx.push_back(index);
			sin >> index;
			flipTetraIdx.push_back(index);
			// get new tetra quality
			sin >> Q_new[0] >> Q_new[1];
			sort(Q_new, Q_new + 2);
			quality_new.push_back(Q_new[0]);
			quality_new.push_back(Q_new[1]);

			text.clear();
			getline(in_quality, text);
		}
	}
	else if (s_flipType == "FLIP23")
	{
		radix_old = 2;
		radix_new = 3;
		Q_new = new double[3];
		getline(in_quality, text);
		while (text != "\n" && text.size() != 0)
		{
			sin.clear();
			sin.str(text);
			// get label
			sin >> label;
			label_quality.push_back(label);
			// get flip tetra index
			sin >> index;
			flipTetraIdx.push_back(index);
			sin >> index;
			flipTetraIdx.push_back(index);
			// get new tetra quality
			sin >> Q_new[0] >> Q_new[1] >> Q_new[2];
			sort(Q_new, Q_new + 3);
			quality_new.push_back(Q_new[0]);
			quality_new.push_back(Q_new[1]);
			quality_new.push_back(Q_new[2]);

			text.clear();
			getline(in_quality, text);
		}
	}
	else if (s_flipType == "FLIP32")
	{
		radix_old = 3;
		radix_new = 2;
		Q_new = new double[2];
		getline(in_quality, text);
		while (text != "\n" && text.size() != 0)
		{
			sin.clear();
			sin.str(text);
			// get label
			sin >> label;
			label_quality.push_back(label);
			// get flip tetra index
			sin >> index;
			flipTetraIdx.push_back(index);
			sin >> index;
			flipTetraIdx.push_back(index);
			sin >> index;
			flipTetraIdx.push_back(index);
			// get new tetra quality
			sin >> Q_new[0] >> Q_new[1];
			sort(Q_new, Q_new + 2);
			quality_new.push_back(Q_new[0]);
			quality_new.push_back(Q_new[1]);

			text.clear();
			getline(in_quality, text);
		}
	}
}

void DataCompare::data_with_different_threshold()
{
	double threshold;
	totalFlipTime = label_quality.size();

	outfile.open("F:\\dataProcessResult\\data_with_threshold", std::ofstream::app);
	threshold = 0.0;
	data_with_threshold(threshold);
	threshold = 0.1;
	data_with_threshold(threshold);
	threshold = 0.2;
	data_with_threshold(threshold);
	threshold = 0.3;
	data_with_threshold(threshold);
	threshold = 0.5;
	data_with_threshold(threshold);
	outfile.close();
}

void DataCompare::data_with_threshold(double threshold)
{
	std::vector<int> label_predict_T;  // 加阈值后的label集合
	int num_out = 0, num_in = 0;  // num_out: 在阈值范围外的记录次数  num_in: 在阈值范围内的记录次数
	int num_do = 0, num_not = 0;  // num_do: 在阈值范围内翻转的次数  num_not：在阈值范围内不翻转的次数
	// classification with threshold
	quality_flip.clear();
	quality_svm.clear();
	quality_flip = quality_svm = quality;
	quality_flip_flag = new int[quality.size()];
	quality_svm_flag = new int[quality.size()];
	for (unsigned int i = 0; i < quality.size(); i++)
	{
		quality_flip_flag[i] = 0;
		quality_svm_flag[i] = 0;
	}
	// calculate quality_flip
	for (unsigned int i = 0; i < label_quality.size(); i++)
	{
		if (label_quality[i] == 1)
			updateTetraQuality(i, quality_flip, quality_flip_flag);
	}
	// calculate quality_svm
	for (unsigned int i = 0; i < predict_0.size(); i++)
	{
		if (!(fabs(predict_0[i] - predict_1[i]) < threshold))
		{
			++ num_in;
			if (label_quality[i] == 1)
				++ num_do;
			else
				++ num_not;

			if (label_predict[i] == 1)
			{
				label_predict_T.push_back(1);
				// update quality vector
				updateTetraQuality(i, quality_svm, quality_svm_flag);
			}
			else
			{
				label_predict_T.push_back(0);
			}
		}
		else
		{
			++ num_out;
			if (label_quality[i] == 1)
			{
				label_predict_T.push_back(1);
				// update quality vector
				updateTetraQuality(i, quality_svm, quality_svm_flag);
			}
			else
			{
				label_predict_T.push_back(0);
			}
		}
	}

	accuracyRate = accuracyRateCalculate(label_quality, label_predict_T);
	do_predict_not = rate_do_predict_not(label_quality, label_predict_T);
	not_predict_do = rate_not_predict_do(label_quality, label_predict_T);
	QualityDistrbtStatistic(quality, oldQD);
	QualityDistrbtStatistic(quality_flip, flipQD);
	QualityDistrbtStatistic(quality_svm, svmQD);

	outfile << entityName.c_str() << "\t" << std::setprecision(0) << threshold * 100 << "%" << "\t";
	outfile << totalFlipTime << "\t" << num_out << "\t" << num_in << "\t" << num_do << "\t" << num_not << "\t";
	outfile << do_predict_not * 1.0 /num_do << "(" << do_predict_not << "/" << num_do << ")\t";
	outfile << not_predict_do * 1.0 /num_not << "(" << not_predict_do << "/" << num_not << ")\t";
	outfile << oldQD[0] << "\t" << oldQD[1] << "\t" << oldQD[2] << "\t" << oldQD[3] << "\t" << oldQD[4] << "\t";
	outfile << flipQD[0] << "\t" << flipQD[1] << "\t" << flipQD[2] << "\t" << flipQD[3] << "\t" << flipQD[4] << "\t";
	outfile << svmQD[0] << "\t" << svmQD[1] << "\t" << svmQD[2] << "\t" << svmQD[3] << "\t" << svmQD[4] << "\n";

	delete [] quality_flip_flag;
	delete [] quality_svm_flag;
}

// 质量分布统计
void DataCompare::QualityDistrbtStatistic(std::vector<double> qualityArray, int * arrayQD)
{
	for (int i = 0; i < 5; i ++)
	{
		arrayQD[i] = 0;
	}

	for (unsigned int i = 0; i < qualityArray.size(); i++)
	{
		if (qualityArray[i] >= 0 && qualityArray[i] < 0.1)
		{
			++ arrayQD[0];
		}
		else if (qualityArray[i] >= 0.1 && qualityArray[i] < 0.3)
		{
			++ arrayQD[1];
		}
		else if (qualityArray[i] >= 0.3 && qualityArray[i] < 0.5)
		{
			++ arrayQD[2];
		}
		else if (qualityArray[i] >= 0.5 && qualityArray[i] < 0.7)
		{
			++ arrayQD[3];
		}
		else if (qualityArray[i] >= 0.7 && qualityArray[i] <= 1.0)
		{
			++ arrayQD[4];
		}
	}
}

double DataCompare::accuracyRateCalculate(std::vector<int> standardVec, std::vector<int> testVec)
{
	if (standardVec.size() != testVec.size())
		return -1;

	int correctTime = 0;
	for (unsigned int i = 0; i < standardVec.size(); i++)
	{
		if (standardVec[i] == testVec[i])
			++ correctTime;
	}
	return correctTime*1.0/standardVec.size();
}

//应翻转，但svm预测为不翻转的情况统计
int DataCompare::rate_do_predict_not(std::vector<int> standardVec, std::vector<int> testVec)
{
	if (standardVec.size() != testVec.size())
		return -1;

	int count = 0;
	for (unsigned int i = 0; i < standardVec.size(); i++)
	{
		if (standardVec[i] == 1 && testVec[i] == 0)
			++ count;
	}
	return count;
}

//不应翻转，但svm预测为翻转的情况统计
int DataCompare::rate_not_predict_do(std::vector<int> standardVec, std::vector<int> testVec)
{
	if (standardVec.size() != testVec.size())
		return -1;

	int count = 0;
	for (unsigned int i = 0; i < standardVec.size(); i++)
	{
		if (standardVec[i] == 0 && testVec[i] == 1)
			++ count;
	}
	return count;
}

// update tetra quality
void DataCompare::updateTetraQuality(int idx, std::vector<double> &qualityV, int *v_flag)
{
	if (flipType == 0x0001)  // flip22
	{
		if (qualityV[flipTetraIdx[idx * radix_old]] < qualityV[flipTetraIdx[idx * radix_old + 1]])
		{
			if (v_flag[flipTetraIdx[idx*radix_old]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old]] == 1 && qualityV[flipTetraIdx[idx*radix_old]] < quality_new[idx*radix_new]))
			{
				qualityV[flipTetraIdx[idx * radix_old]] = quality_new[idx * radix_new];
				v_flag[flipTetraIdx[idx*radix_old]] = 1;
			}
			if (v_flag[flipTetraIdx[idx*radix_old+1]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old+1]] == 1 && qualityV[flipTetraIdx[idx*radix_old+1]] < quality_new[idx*radix_new+1]))
			{
				qualityV[flipTetraIdx[idx * radix_old + 1]] = quality_new[idx * radix_new + 1];
				v_flag[flipTetraIdx[idx*radix_old + 1]] = 1;
			}
		}
		else
		{
			if (v_flag[flipTetraIdx[idx*radix_old]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old]] == 1 && qualityV[flipTetraIdx[idx*radix_old]] < quality_new[idx*radix_new+1]))
			{
				qualityV[flipTetraIdx[idx * radix_old]] = quality_new[idx * radix_new + 1];
				v_flag[flipTetraIdx[idx*radix_old]] = 1;
			}
			if (v_flag[flipTetraIdx[idx*radix_old+1]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old+1]] == 1 && qualityV[flipTetraIdx[idx*radix_old+1]] < quality_new[idx*radix_new]))
			{
				qualityV[flipTetraIdx[idx * radix_old + 1]] = quality_new[idx * radix_new];
				v_flag[flipTetraIdx[idx*radix_old + 1]] = 1;
			}
			
		}
	}
	else if (flipType == 0x0002)  // flip23
	{
		if (qualityV[flipTetraIdx[idx * radix_old]] < qualityV[flipTetraIdx[idx * radix_old + 1]])
		{
			if (v_flag[flipTetraIdx[idx*radix_old]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old]] == 1 && qualityV[flipTetraIdx[idx*radix_old]] < quality_new[idx*radix_new]))
			{
				qualityV[flipTetraIdx[idx * radix_old]] = quality_new[idx * radix_new];
				v_flag[flipTetraIdx[idx*radix_old]] = 1;
			}
			if (v_flag[flipTetraIdx[idx*radix_old+1]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old+1]] == 1 && qualityV[flipTetraIdx[idx*radix_old+1]] < quality_new[idx*radix_new+1]))
			{
				qualityV[flipTetraIdx[idx * radix_old + 1]] = quality_new[idx * radix_new + 1];
				v_flag[flipTetraIdx[idx*radix_old + 1]] = 1;
			}
			qualityV.push_back(quality_new[idx * radix_new + 2]);
		}
		else
		{
			if (v_flag[flipTetraIdx[idx*radix_old]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old]] == 1 && qualityV[flipTetraIdx[idx * radix_old]] < quality_new[idx * radix_new + 1]))
			{
				qualityV[flipTetraIdx[idx * radix_old]] = quality_new[idx * radix_new + 1];
				v_flag[flipTetraIdx[idx*radix_old]] = 1;
			}
			if (v_flag[flipTetraIdx[idx*radix_old+1]] == 0 ||
				(v_flag[flipTetraIdx[idx*radix_old+1]] == 1 && qualityV[flipTetraIdx[idx * radix_old + 1]] < quality_new[idx * radix_new]))
			{
				qualityV[flipTetraIdx[idx * radix_old + 1]] = quality_new[idx * radix_new];
				v_flag[flipTetraIdx[idx*radix_old + 1]] = 1;
			}
			qualityV.push_back(quality_new[idx * radix_new + 2]);
		}
	}
	else if (flipType == 0x0004)  // flip32
	{
		qualityV[flipTetraIdx[idx * radix_old]] = -1;
		qualityV[flipTetraIdx[idx * radix_old + 1]] = -1;
		qualityV[flipTetraIdx[idx * radix_old + 2]] = -1;
		qualityV.push_back(quality_new[idx * radix_new]);
		qualityV.push_back(quality_new[idx * radix_new + 1]);
	}
}