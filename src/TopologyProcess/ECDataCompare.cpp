#include <TopologyProcess/ECDataCompare.h>
#include <Windows.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>

void ECDataCompare::ec_data_compare()
{
	if (!initialize())
		return;

	if (is_threshold)
	{
		data_analyse_threshold(0.00);
		data_analyse_threshold(0.05);
		data_analyse_threshold(0.10);
		data_analyse_threshold(0.15);
		data_analyse_threshold(0.20);
		data_analyse_threshold(0.25);
		data_analyse_threshold(0.30);
	}
	else
		data_analyse();
}

void ECDataCompare::data_analyse()
{
	do_predict_not = not_predict_do = do_predict_do = not_predict_not = 0;
	for (int i = 0; i < ecTime; i++)
	{
		if (ECFlag[i] == 1 && predictVec[i] == 0)
			++ do_predict_not;
		if (ECFlag[i] == 0 && predictVec[i] == 1)
			++ not_predict_do;
		if (ECFlag[i] == 1 && predictVec[i] == 1)
			++ do_predict_do;
		if (ECFlag[i] == 0 && predictVec[i] == 0)
			++ not_predict_not;
	}

	//QImproveQSDetail();
	std::string fn = path;
	fn += "ec_compare_result";
	outfile.open(fn.c_str(), std::ofstream::app);
	if (!outfile)
		return;
	outfile << modelName.c_str() << "\t" << ecTime << "\t" << succ_time << "\t" 
		<< std::setprecision(6) << (do_predict_do+not_predict_not)*1.0/ecTime*100 << "%" << "\t"
		<< do_predict_do << "(" << do_predict_do * 1.0/succ_time << ")" << "\t" 
		<< not_predict_not << "(" << not_predict_not * 1.0/(ecTime-succ_time) << ")" << "\t" 
		<< do_predict_not << "(" << do_predict_not * 1.0/succ_time << ")" << "\t" 
		<< not_predict_do << "(" << not_predict_do * 1.0/(ecTime-succ_time) << ")" << "\t";
	outfile << "\n";
	outfile.close();
}

void ECDataCompare::data_analyse_threshold(double threshold)
{
	int predictCount, unpredictCount;
	predictCount = unpredictCount = 0;
	do_predict_not = not_predict_do = do_predict_do = not_predict_not = 0;
	for (int i = 0; i < ecTime; i++)
	{
		if (ECFlag[i] == 1)
		{
			if ( fabs(predictRate[i*2+1]-predictRate[i*2]) >= threshold )
			{
				++ predictCount;
				if (predictVec[i] == 1)
					++do_predict_do;
				else
					++do_predict_not;
			}
			else
				++ unpredictCount;
		}
		if (ECFlag[i] == 0)
		{
			if ( fabs(predictRate[i*2+1]-predictRate[i*2]) >= threshold )
			{
				++ predictCount;
				if (predictVec[i] == 0)
					++ not_predict_not;
				else
					++ not_predict_do;
			}
			else
				++ unpredictCount;
		}
	}

	std::string fn = path;
	fn += "ec_compare_result_threshold";
	outfile.open(fn.c_str(), std::ofstream::app);
	if (!outfile)
		return;
	outfile << modelName.c_str() << "\t" << ecTime << "\t" << succ_time << "\t" 
		<< std::setprecision(0) << threshold * 100 << "%" << "\t" 
		<< predictCount << "(" << std::setprecision(6) << predictCount*1.0/ecTime*100 << "%)" << "\t"
		<< unpredictCount << "(" << std::setprecision(6) << unpredictCount*1.0/ecTime*100 << "%)" << "\t"
		<< std::setprecision(6) << (do_predict_do+not_predict_not)*1.0/predictCount*100 << "%" << "\t"
		<< do_predict_do << "(" << std::setprecision(6) << do_predict_do*1.0/ecTime*100 << "%)" << "\t" 
		<< not_predict_not << "(" << std::setprecision(6) << not_predict_not*1.0/ecTime*100 << "%)" << "\t" 
		<< do_predict_not << "(" << std::setprecision(6) << do_predict_not*1.0/ecTime*100 << "%)" << "\t" 
		<< not_predict_do << "(" << std::setprecision(6) << not_predict_do*1.0/ecTime*100 << "%)" << "\t";

	outfile << "\n";
	outfile.close();
}

int ECDataCompare::initialize()
{
	// get model name
	std::string fileName = ecfile;
	char *tokenPtr = strtok((char *)(fileName.c_str()), "\/, \\");
	while(tokenPtr != NULL)
	{
		fileName = std::string(tokenPtr);
		tokenPtr = strtok(NULL, "\/, \\");
		if (tokenPtr != NULL)
		{
			path += fileName;
			path += "\\";
		}
	}
	tokenPtr = strtok((char*)(fileName.c_str()), "_");  // ec_
	tokenPtr = strtok(NULL, "_");
	modelName = std::string(tokenPtr);
	tokenPtr = strtok(NULL, "_");
	while(tokenPtr != NULL && std::string(tokenPtr) != "svm")
	{
		modelName += "_";
		modelName += std::string(tokenPtr);
		tokenPtr = strtok(NULL, "_");
	}

	std::ifstream infile;
	infile.open(ecfile.c_str());
	if (!infile)
	{
		MessageBox(0,TEXT("EC file is not exist!"),TEXT("Error"),MB_OK);
		return 0;
	}
	infile.close();
	infile.open(pfile.c_str());
	if (!infile)
	{
		MessageBox(0,TEXT("Predict file is not exist!"),TEXT("Error"),MB_OK);
		return 0;
	}
	infile.close();

	// can not change the file reading order
	read_svm_file();
	read_predict_file();
	return 1;
}
void ECDataCompare::read_predict_file()
{
	std::ifstream infile;
	infile.open(pfile.c_str());

	std::string str;
	int intData;
	double doubleData;
	if (is_threshold)
	{
		int idx = 0;
		predictVec.resize(ecTime, 0);
		predictRate.resize(ecTime*2, 0.0);

		//获取0和1类的百分比的顺序
		getline(infile, str);
		if (str != "")
		{
			int t1,t2;
			sscanf(str.c_str(), "%*s %d %d", &t1, &t2);
			if (t1 == 0)
				while (getline(infile, str) && str!="")
				{
					sscanf(str.c_str(),"%d %lf %lf", &predictVec[idx], &predictRate[idx*2], &predictRate[idx*2+1]);
					++idx;
				}
			else if (t1 == 1)
				while (getline(infile, str) && str!="")
				{
					sscanf(str.c_str(),"%d %lf %lf", &predictVec[idx], &predictRate[idx*2+1], &predictRate[idx*2]);
					++idx;
				}
		}
	}
	else
	{
		int idx = 0;
		predictVec.resize(ecTime, 0);
		while (getline(infile, str) && str!="")
		{
			if (str[0] != '0' && str[1] != '1')
				continue;
			sscanf(str.c_str(),"%d%*s", &predictVec[idx++]);
		}
	}

	infile.close();
}

//不能有特征缺失
void ECDataCompare::read_svm_file()
{
	std::ifstream infile;
	infile.open(ecfile.c_str());

	std::string str;
	getline(infile, str);
	if (str.size() < 2)
		return;

	int iData;
	double dData;
	charaNum = 0;
	succ_time = 0;
	ECFlag.clear();
	ECCharacter.clear();

	//read the first line of svm file
	char *tokenPtr = strtok((char *)(str.c_str()), " ");
	//get ec flag
	sscanf(tokenPtr,"%d",&iData);
	ECFlag.push_back(iData);
	if (iData)
		++succ_time;
	//get ec character
	tokenPtr = strtok(NULL, " ");
	while(tokenPtr != NULL)
	{
		++ charaNum;
		sscanf(tokenPtr, "%d:%lf", &iData, &dData);
		ECCharacter.push_back(dData);
		tokenPtr = strtok(NULL, " ");
	}

	//read other lines
	while(getline(infile, str) && str.size() > 2)
	{
		char *tokenPtr = strtok((char *)(str.c_str()), " ");
		//get ec flag
		sscanf(tokenPtr,"%d",&iData);
		ECFlag.push_back(iData);
		if (iData)
			++succ_time;
		//get ec character
		tokenPtr = strtok(NULL, " ");
		while(tokenPtr != NULL)
		{
			sscanf(tokenPtr, "%d:%lf", &iData, &dData);
			ECCharacter.push_back(dData);
			tokenPtr = strtok(NULL, " ");
		}
	}
	infile.close();
	ecTime = ECFlag.size();
}