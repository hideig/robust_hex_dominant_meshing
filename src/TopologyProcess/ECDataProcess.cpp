#include <TopologyProcess/ECDataProcess.h>
#include <stdlib.h> 
#include <string>
#include <sstream>

void ECDataProcess::initialize()
{
	std::ifstream infile;
	infile.open(fileName.c_str());
	if (!infile)
		return;

	std::string str;
	char *tokenPtr = NULL;
	//数据总数，成功次数
	getline(infile, str);
	//infile >> line_counter >> succ_num;
	sscanf(str.c_str(), "%d %d", &line_counter, &succ_num);

	ECFlag.resize(line_counter, 0);
	ETVec.resize(line_counter*3, 0);
	QVec.resize(line_counter*4, 0);
	locationVec.resize(line_counter, 0);
	QImproveVec.resize(line_counter, 0);
	sideRatio.resize(line_counter*9, 0);
	minSineVec.resize(line_counter, 0);
	squareRoot.resize(line_counter, 0);
	jacValAroudEdge.resize(line_counter);

	int intData;
	double doubleData;
	for (int i = 0; i < line_counter; i++)
	{
		getline(infile, str);

		//操作标志
		tokenPtr = strtok((char *)(str.c_str()), "\t");
		sscanf(tokenPtr,"%d",&intData);
		//infile >> intData;
		ECFlag[i] = intData;

		//两端点周围网格数
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%d",&intData);
		//infile >> intData;
		ETVec[i*3] = intData;

		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%d",&intData);
		//infile >> intData;
		ETVec[i*3+1] = intData;

		//边周围网格数
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%d",&intData);
		//infile >> intData;
		ETVec[i*3+2] = intData;

		//总网格数
		//infile >> intData;
		tokenPtr = strtok(NULL, "\t");

		// 操作前网格的最好、最差质量值
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		QVec[i*4] = doubleData;

		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		QVec[i*4+1] = doubleData;

		//操作后的网格数
		//infile >> intData;
		tokenPtr = strtok(NULL, "\t");

		//操作后网格的最好、最差质量值
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		QVec[i*4+2] = doubleData;

		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		QVec[i*4+3] = doubleData;

		//边收缩前最差质量网格所处的位置（0：边周围  1：端点周围）
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%d",&intData);
		//infile >> intData;
		locationVec[i] = intData;

		//最差网格质量提高的绝对值
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		QImproveVec[i] = doubleData;

		//最长最短边之比
		for (int k = 0; k < 9; k ++)
		{
			tokenPtr = strtok(NULL, "\t");
			sscanf(tokenPtr,"%lf",&doubleData);
			//infile >> doubleData;
			sideRatio[i*9+k] = doubleData;
		}

		//最小二面角sine值
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		minSineVec[i] = doubleData;

		//square root 质量值
		tokenPtr = strtok(NULL, "\t");
		sscanf(tokenPtr,"%lf",&doubleData);
		//infile >> doubleData;
		squareRoot[i] = doubleData;

		//边周围网格顶点的雅克比值
		std::vector<double> tVec;
		tokenPtr = strtok(NULL, "\t");
		while(tokenPtr != NULL)
		{
			sscanf(tokenPtr, "%lf", &doubleData);
			tVec.push_back(doubleData);
			tokenPtr = strtok(NULL, "\t");
		}
		jacValAroudEdge[i] = tVec;

	}
	infile.close();

	// initialize array
	for (int i = 0; i < 15; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			TetraRange_QIWT[i][j] = 0;
		}
	}

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 18; j++)
		{
			TetraRange_QIET[i][j] = 0;
		}
	}
}

void ECDataProcess::edge_contract_data_process()
{
	initialize();
	save_svm_data();
	//tetra_range_QIWT();
	//tetra_range_QIET();
}

void ECDataProcess::tetra_range_QIWT()
{
	int i, j;
	for (int iter = 0; iter < line_counter; iter++)
	{
		array_idx_TetraRange_QIWT(QImproveVec[iter], QVec[iter*4+1], i, j);
		if ((i>-1 && i<15) && (j>-1 && j<10))
			++ TetraRange_QIWT[i][j];
	}

	std::ofstream outfile;
	std::string fn = fileName;
	fn += "_TetraRange_QIWT";
	outfile.open(fn.c_str());
	int rowTotal;
	for (i = 0; i < 15; i++)
	{
		rowTotal = 0;
		for (j = 0; j < 10; j++)
		{
			rowTotal += TetraRange_QIWT[i][j];
			outfile << TetraRange_QIWT[i][j] << "\t";
		}
		outfile << rowTotal << "\n";
	}
	outfile.close();
}

void ECDataProcess::tetra_range_QIET()
{
	int i, j;
	for (int iter = 0; iter < line_counter; iter++)
	{
		i = j = -1;
		array_idx_TetraRange_QIET(QImproveVec[iter], ETVec[iter*3+2], i, j);
		++ TetraRange_QIET[i][j];
	}
	std::ofstream outfile;
	std::string fn = fileName;
	fn += "_TetraRange_QIET";
	outfile.open(fn.c_str());
	for (i = 0; i < 18; i++)
	{
		for (j = 0; j < 9; j++)
			outfile << TetraRange_QIET[j][i] << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void ECDataProcess::array_idx_TetraRange_QIWT(double row, double column, int &i, int &j)
{
	if (row <= 0)
		i = -1;
	else if (row <= 0.01)
		i = 0;
	else if (row <= 0.02)
		i = 1;
	else if (row <= 0.03)
		i = 2;
	else if (row <= 0.04)
		i = 3;
	else if (row <= 0.05)
		i = 4;
	else if (row <= 0.06)
		i = 5;
	else if (row <= 0.07)
		i = 6;
	else if (row <= 0.08)
		i = 7;
	else if (row <= 0.09)
		i = 8;
	else if (row <= 0.1)
		i = 9;
	else if (row <= 0.2)
		i = 10;
	else if (row <= 0.3)
		i = 11;
	else if (row <= 0.4)
		i = 12;
	else if (row <= 0.5)
		i = 13;
	else if (row <= 0.6)
		i = 14;

	if (column <= 0)
		j = -1;
	else if (column <= 0.1)
		j = 0;
	else if (column <= 0.2)
		j = 1;
	else if (column <= 0.3)
		j = 2;
	else if (column <= 0.4)
		j = 3;
	else if (column <= 0.5)
		j = 4;
	else if (column <= 0.6)
		j = 5;
	else if (column <= 0.7)
		j = 6;
	else if (column <= 0.8)
		j = 7;
	else if (column <= 0.9)
		j = 8;
	else if (column <= 1.0)
		j = 9;
}


void ECDataProcess::array_idx_TetraRange_QIET(double row, int column, int &i, int &j)
{
	if (row <= 0)
		i = -1;
	else if (row <= 0.03)
		i = 0;
	else if (row <= 0.05)
		i = 1;
	else if (row <= 0.07)
		i = 2;
	else if (row <= 0.1)
		i = 3;
	else if (row <= 0.2)
		i = 4;
	else if (row <= 0.3)
		i = 5;
	else if (row <= 0.4)
		i = 6;
	else if (row <= 0.5)
		i = 7;
	else if (row <= 0.6)
		i = 8;
	else
		i = -1;

	if (column >= 3 && column <= 20)
		j = column - 3;
	else
		j = -1;
}

void ECDataProcess::save_svm_data()
{
	std::string file;
	std::string file_sort_TN[9];
	file = fileName;

	if (threshold != 0)
	{
		char *ts;
		int dec, sign;
		int ndig = 2;
		ts = ecvt(threshold, ndig, &dec, &sign); 
		file += "_";
		file += std::string(ts);
	}

	if (svm_tetra_num != 0x0000)
	{
		char tchar[5];
		for (int i = 0; i < 9; i ++)
		{
			file_sort_TN[i] = file;
			itoa(i+3, tchar, 10);
			file_sort_TN[i] += "_TN";
			file_sort_TN[i] += std::string(tchar);
			file_sort_TN[i] += "_svm";
			ofSortTN[i].open(file_sort_TN[i].c_str());
			if (!ofSortTN[i])
				return;
		}
	}
	else
	{
		file += "_svm";
		outfile.open(file.c_str());
		if (!outfile)
			return;
	}

	int ec_flag;
	std::vector<double> attr;
	for (int i = 0; i < line_counter; i++)
	{
		attr.clear();
		// 边周围网格的最差质量值
		if (svm_attr & 0x0001)
			attr.push_back(QVec[i*4+1]);
		// square root
		if (svm_attr & 0x0020)
			attr.push_back(squareRoot[i]);
		// 边周围网格最差二面角(共享收缩边的两个面的二面角)
		if (svm_attr & 0x0002)
			attr.push_back(minSineVec[i]);
		// 边周围网格数
		if (svm_attr & 0x0004)
			attr.push_back(3.0/ETVec[i*3+2]);
		// 边周围网格的长短边比
		if (svm_attr & 0x0008)
			attr.push_back(sideRatio[i*9+4]);
		// 端点周围网格的长短边比
		if (svm_attr & 0x0010)
			attr.push_back(sideRatio[i*9+7]);

		// 收缩边周围的网格端点的雅克比值
		if (svm_tetra_num != 0x0000)
		{
			if (ETVec[i*3+2]<11)
				for (unsigned int k = 0; k < jacValAroudEdge[i].size(); k++ )
					attr.push_back(jacValAroudEdge[i][k]);
			else
				for (unsigned int k = 0; k < 22; k++)
					attr.push_back(jacValAroudEdge[i][k]);
		}

		if (QImproveVec[i] <= threshold)
			ec_flag = 0;
		else
			ec_flag = 1;

		if (svm_tetra_num != 0x0000)
		{
			if (ETVec[i*3+2] > 2 && ETVec[i*3+2] < 11)
				save_svm_data_single_line(ec_flag, attr, ofSortTN[ETVec[i*3+2]-3]);
			else if (ETVec[i*3+2] >= 11)
				save_svm_data_single_line(ec_flag, attr, ofSortTN[8]);
		}
		else
			save_svm_data_single_line(ec_flag, attr, outfile);
	}
	outfile.close();
	if (svm_tetra_num != 0x0000)
	{
		for (int i = 0; i < 9; i ++)
		{
			ofSortTN[i].close();
		}
	}
}

void ECDataProcess::save_svm_data_single_line(int succ, std::vector<double> attr, std::ofstream &outfile)
{
	// flag
	outfile << succ << " ";
	// attributes
	for (unsigned int i = 0; i < attr.size(); i++)
	{
		outfile << i+1 << ":" << attr[i] << " ";
	}
	outfile << "\n";
}