#ifndef _TETRA_MESH_TOOL_JOURNAL_
#define _TETRA_MESH_TOOL_JOURNAL_
#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/top.h>
#include <deque>
namespace VolumeMesh
{
	/* a struct for journal record */
	struct Record
	{
		int type;
		int optclass;
		std::vector<Point> pointvec;
		std::vector<TetraHandle> thvec;
		std::vector<PointHandle> phvec;
		std::vector<VertexHandle> vhvec;
	};

	typedef std::deque<Record> JournalType;

	class Journal
	{
	public:
		Journal(){}
		~Journal(){}

		bool empty()
		{
			return journallist.size() == 0;
		}

		void clear()
		{
			journallist.clear();
		}

		std::deque<Record>::size_type size()
		{
			return journallist.size();
		}

		Record topJournal()
		{
			if (journallist.size() == 0)
				return Record();
			return journallist.back();
		}

		Record popJournal()
		{
			if (journallist.size() == 0)
				return Record();
			Record rc = journallist.back();
			journallist.pop_back();
			return rc;
		}

		/* insert a record to a journal container */
		void pushJournal(int type, const std::vector<Point> &pointvec, const std::vector<TetraHandle> &thvec_,  
			             const std::vector<PointHandle> &phvec, const std::vector<VertexHandle> &vhvec_);


	protected:
	private:
		std::deque<Record> journallist;
	};
}
#endif