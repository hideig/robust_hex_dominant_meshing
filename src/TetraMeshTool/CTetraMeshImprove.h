#ifndef _CUDA_TETRA_MESH_IMPROVE_
#define _CUDA_TETRA_MESH_IMPROVE_

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <TetraMeshTool/top.h>
#include <TetraMeshTool/QualityCalculator.h>
#include <map>
#include <fstream>
using namespace std;

namespace VolumeMesh
{
	struct cu_face;
	// cuda interface
	extern "C" void cuda_vertexSmoothing(float *points, int pointcnt, int *hneighbour, int *neighbourcnt, int *pointgroup, int pointgroupcnt, 
		float *newpoints, int *hincidenttet, int *incidenttetcnt, int *meshtets, int tetcnt, 
		int largestn, int largesttet, int smoothpasscnt, float &time);

	extern "C" void cuda_tetquality(float *points, int pointcnt, int *meshtets, int tetcnt, int qualmeasure, float &minqual);

	extern "C" void cuda_flip23(float *points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
		                        int *halfface, int halffacecnt, int qualmeasure, float& qualbefore, float& qualafter, 
								int &flipsucc, float &time);

	extern "C" void cuda_flip32(float *points, int pointcnt, int *meshtets_, int tetcnt, int *flipedge, int edgecnt, 
		                        int qualmeasure, float& qualbefore_, float& qualafter_, int &flipsucc, float &time);

	extern "C" void cuda_flip32_new(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
	                                int *halfedge, int halfedgecnt, int qualmeasure, float &qualbefore_, float &qualafeter_,
		                            int &flipsucc, float &time);

	extern "C" void cuda_newflip(float *points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
		                         int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
		                         int &flipsucc, float &time);

	extern "C" void cuda_edgeContraction(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
		                                 int *halfedge, int halfedgecnt, int *incidenttet, int* incidenttetcnt, int &largesttetcnt, 
										 int qualmeasure, float &qualbefore_, float &qualafeter_, int &succ, float &time);

	extern "C" void cuda_vertexInsertion(float *&points, int pointcnt, int *meshtets_, int tetcnt, int *face, int facecnt, 
		                                 int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
		                                 int &succ, float &time);


	struct cu_halfface;
	struct cu_face
	{
		int hf[2];
		float quality;
		float val;
	};

	struct cu_halfface
	{
		int face;
		int pointhandle[3];
	};

	struct cu_edge
	{
		bool is_boundary;
		int v[2];
		int halfedgecnt;
		int *halfedge;
		float quality;
		float val;
	};

	struct cu_halfedge
	{
		int edgehandle;
		int fromv;
		int tov;
	};

	struct cu_tetra
	{
		bool isboundary;
		int v[4];
		float quality;
		float val;
		int fliptype;
		int newflipvertex;
		int newflipface[4];
	};

	template<class T>
	struct queue
	{
		queue()
		{
			initialize();
		}

		~queue()
		{
			node *tmp;
			tmp=head;
			while(tmp!=NULL)
			{
				head=head->next;
				delete tmp;
				tmp=head;
			}
		}

		void initialize()
		{
			head=new node(-1);
			cur=head;
			len=0;
		}

		bool empty() const
		{
			return len==0;
		}

		T& back()
		{
			return cur->val;
		}

		const T& back() const
		{
			return back();
		}

		void pop()
		{
			if(head->next==cur)
			{
				delete head->next;
				head->next=NULL;
			}else
			{
				node* tmp=head->next;
				head->next=tmp->next;
				delete tmp;
			}
			--len;
		}

		T& front()
		{
			return head->next->val;
		}

		const T& front() const
		{
			return front();
		}

		void push(const T& val)
		{
			node *tmp=new node(val);
			if(!len)
				head->next = tmp;
			cur->next=tmp;
			cur=tmp;
			++len;
		}

		int size()
		{
			return len;
		}

		typedef struct node1
		{
			node1 *next;
			T val;
			node1(T v):val(v),next(NULL){}
		}node;

		int len;
		node *head;
		node *cur;

	};

	class CTetraMeshImprove
	{
	public:
		CTetraMeshImprove():mesh(NULL){}
		CTetraMeshImprove(TetraMesh *mesh_):mesh(mesh_){}
		~CTetraMeshImprove(){mesh = NULL;}

		void setMesh(TetraMesh *mesh_)
		{
			mesh = mesh_;
		}

		void pointsandneighbours(float* points, int pointcnt, int** &neighbours, int* neighbourcnt, 
			                     int** &incidenttet, int* incidenttetcnt, float &time);
		void pointsNeighbourTetra(TetraMesh *mesh, int** &incidenttet, int* &incidenttetcnt, int &largesttetcnt, float &time);

		void getcudapoints(TetraMesh *mesh, float* points);

		void cuda_edgeContraction_test(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
			int *halfedge, int halfedgecnt, int *incidenttet, int* incidenttetcnt, int &largesttetcnt, 
			int qualmeasure, float &qualbefore_, float &qualafeter_, int &succ, float &time);
// 
// 		void edge_contraction_test(float *points, int *meshtets, cu_edge *edge, cu_halfedge *halfedge,
// 			int *selectedge, int selectedgecnt);
// 		void edge_contraction_explore_test(float *points, int *meshtets, cu_edge *edge, int edgecnt, cu_halfedge *halfedge, 
// 			int *incidenttet, int *incidenttetcnt, int largesttetcnt, int qualmeasure);
// 		void getECNewTetras_test(int *newtetra, int *newtetracnt, int *fromtetra, int fromtetracnt, int *totetra, int totetracnt,
// 			int *edgestar, int edgestarcnt);

	protected:
	private:
		TetraMesh *mesh;
	};


}
#endif