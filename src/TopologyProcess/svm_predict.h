#ifndef TOPOLOGY_PROCESS_SVM_PREDICT
#define TOPOLOGY_PROCESS_SVM_PREDICT
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "svm.h"
class SVMPredict
{
public:
	SVMPredict()
	{
		max_nr_attr = 64;
		predict_probability=0;
		line = NULL;
	}
	~SVMPredict(){}
protected:
private:
	struct svm_node *x;
	int max_nr_attr;

	struct svm_model* model;
	int predict_probability;

	static char *line;
	static int max_line_len;
};

#endif