#include <TetraMeshTool/Journal.h>

namespace VolumeMesh
{
	/* insert a record to a journal container */
	void Journal::pushJournal(int type, const std::vector<Point> &pointvec_, const std::vector<TetraHandle> &thvec_,  
		                      const std::vector<PointHandle> &phvec_, const std::vector<VertexHandle> &vhvec_)
	{
		Record record;
		/* copy information */
		record.type = type;
		copy(pointvec_.begin(), pointvec_.end(), back_inserter(record.pointvec));
		copy(thvec_.begin(), thvec_.end(), back_inserter(record.thvec));
		copy(phvec_.begin(), phvec_.end(), back_inserter(record.phvec));
		copy(vhvec_.begin(), vhvec_.end(), back_inserter(record.vhvec));

		/* type specific stuff */
		switch (type)
		{
		case SMOOTHVERTEX:
			record.optclass = SMOOTH;
			break;
		case INSERTTET:
			record.optclass = TOPOLOGICAL;
			break;
		case DELETETET:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP23:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP22:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP32:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP13:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP12:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP14:
			record.optclass = TOPOLOGICAL;
			break;
		case FLIP41:
			record.optclass = TOPOLOGICAL;
			break;
		case EDGECONTRACT:
			record.optclass = TOPOLOGICAL;
			break;
		case INSERTVERTEX:
			record.optclass = INSERTDELETE;
			break;
		case CLASSIFY:
			record.optclass = LABEL;
			break;
		default:
			printf("Undefined operation type %d\n", record.type);
			exit(1);
		};

		journallist.push_back(record);
	}

}