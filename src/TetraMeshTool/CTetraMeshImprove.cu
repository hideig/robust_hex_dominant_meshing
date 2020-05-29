
#include <stdio.h>
#include <stdlib.h> 
#include <cuda.h>   
#include <cuPrintf.cu>
#include <sm_12_atomic_functions.h>


//#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.26
#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.65
//#define VERTEX_INSERTION_QUALITY_THRESHOLD 0.45 
#define MINIMPROVEMENT 1.0e-6

#define MAXCAVITYFACES 30
#define MAXCAVITYTETS 50
#define MAXCAVITYEDGES 150
#define MAXSTACKTET 60
#define MAXSTACKFACE 90

#define HUGEFLOAT  1.0e10
#define HUGEQUAL 1.0e10
#define MINFACING 1.0e-7
#define MINFLIPIMPROVE 1.0e-6
#define CAVDEPTHLIMIT 6
/*vertex insert*/
#define EDGELABEL 0
//#define TETLABEL 1
#define CAVLABEL 0
#define ANTICAVLABEL 1
#define NOLABEL 2
#define DEPTHTABLESIZE 10
#define NOCAVITYTET -1
#define NOCAVITYFACE -1
#define GHOSTTET -1
#define imin(a,b) (a<b?a:b)

const int BlockPerGrid = 4096;
const int ThreadPerBlock = 128;
//const int ThreadsPerBlock = 256;

/************************************************************************/
/*                     cuPrint initialization                           */
/************************************************************************/

bool InitGPUSet()  
{  
	char GPU[100] = "GPU: ";  
	cudaDeviceProp tCard;  
	int num = 0;  
	if(cudaSuccess == cudaGetDeviceCount(&num))  
	{  
		for(int i = 0; i < num; ++ i)  
		{  
			cudaSetDevice(i);  
			cudaGetDeviceProperties(&tCard, i);  
			puts(strcat(GPU , tCard.name));//返回的就是链接后的结果,也为其的嵌套使用提供了条件   
		}  
	}  
	else  return false;  
	return true;  
}  
bool cuPrintInit()  
{  
	cudaError_t err = cudaPrintfInit();  
	if(0 != strcmp("no error", cudaGetErrorString(err)))  return false;  
	return true;  
}  
__global__ void displayGPU_demo()  
{  
	int bsize = blockDim.x;  
	int bid = blockIdx.x;  
	int tid = bid * bsize + threadIdx.x;  
	cuPrintf("当前执行kernel的 block 编号:\t%d\n", bid);  
	cuPrintf("当前执行kernel的 thread 在当前块中编号:\t%d\n", threadIdx.x);  
	cuPrintf("当前执行kernel的 thread 全局编号:\t%d\n", tid);  
	cuPrintf("thread over\n\n");  
}

// 示例：将程序放在最后一个else里即可
extern "C" void testCudaPrintf()
{
	if(!InitGPUSet())  puts("device is not ready!");  
	else if(!cuPrintInit())  puts("device is not ready!");  
	else  
	{  
		displayGPU_demo<<<2, 3>>>();  
		cudaPrintfDisplay(stdout, true);//true输出是哪一个block的第几个thread在执行本条输出语句，形如：[blockID, threadID]；false不输出   
		cudaPrintfEnd();  
	} 
}

/*************************end of cuPrintf initialization ****************************/

/************************************************************************/
/* Data structure                                                       */
/************************************************************************/
struct cu_point
{
	float vec3f[3];
};

struct cu_tet
{
	int v[4];
};


struct cu_halfface
{
	int pointhandle[3];
	int face;
};

struct cu_halfedge
{
	int edgehandle;
	int fromv;
	int tov;
};

struct cu_face
{
	int hf[2];
};

struct cu_flip23face
{
	int hf[2];
	float quality;  // worst quality of incident tets
	float val;      // the improvement of flip23
};



struct CavityFace
{
	int handle;
	float quality;
	int child;
	bool inH;
};

struct CavityEdge
{
	float qual;
	int label;
	int parent;
	int child;
	int childnum;
};

struct CavityTet
{
	int handle;
	float quality;
	int depth;
	CavityFace outfaces[4];
	int outfacesize;
	int parents[3];
	int label;
};

struct cu_InsertTet
{
	int v[4];
	int deletetet[MAXCAVITYTETS];
	int deletetetcnt;
	int cavityface[MAXCAVITYFACES];
	int cavityfacecnt;
	float insertpoint[3];
	float quality;
	float val;
};

struct cu_edge
{
	bool is_boundary;
	int halfedgecnt;
	int halfedge[100];
	float quality;
	float val;
};


struct cu_flip32edge
{
	int p[2];
	int tet[3];
	float quality;  // worst quality of tets
	float val;      // the improvement of flip32
	int order;
};

struct cu_tetra
{
	bool isboundary;
	int v[4];
	float quality;
	float val;
	int fliptype;
	int strategy;
	int newflipvertex;
	int newflipface[4];
	int tet[4];
	int flippoint[2];
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
		while(tmp!=NULL)
		{
			tmp=head;
			head=head->next;
			delete tmp;
			tmp=NULL;
		}
		delete cur;
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
/********************************* end of data structure**************************************/

/****************************** Function Statement *********************************/
extern "C" void cuda_tetquality(float *points, int pointcnt, int *meshtets, int tetcnt, int qualmeasure, float &minqual);

/****************************End of Function Statement******************************/

/************************************************************************/
/*  tetrahedron quality calculate                                       */
/************************************************************************/
/* types of quality measures that may be used */
extern enum CudaTetQualityMetrics
{
	CUDA_QUAL_MINSINE,
	CUDA_QUAL_BIASEDMINSINE,
	CUDA_QUAL_RADIUSRATIO,
	CUDA_QUAL_VLRMS3RATIO,
	CUDA_QUAL_MEANSINE,
	CUDA_QUAL_MINSINEANDEDGERATIO,
	CUDA_QUAL_WARPEDMINSINE,
	CUDA_QUAL_MINANGLE,
	CUDA_QUAL_MAXANGLE
};

__device__ void vector_cross(float *a, float *b, float *c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ float vector_dot(float *a, float *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

__device__ void vector_add(float *a, float *b, float *c)
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

__device__ void vector_minus(float *a, float *b, float *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

__device__ float minsine(float point[4][3])
{
	float t[3], u[3], v[3]; /* tet vectors */
	float temp[3];
    float edgelength[3][4]; /* the lengths of each of the edges of the tet */
    float facenormal[4][3]; /* the normals of each face of the tet */
    float dx, dy, dz;       /* intermediate values of edge lengths */
    float facearea2[4];     /* areas of the faces of the tet */
    float pyrvolume;        /* volume of tetrahedron */
    float sine2, minsine2;  /* the sine (squared) of the dihedral angle */
    int i, j, k, l;          /* loop indices */
    
    /* calculate the volume*6 of the tetrahedron */
	vector_minus(point[1], point[0], t);
	vector_minus(point[2], point[0], u);
	vector_minus(point[3], point[0], v);
	vector_cross(t, u, temp);
    pyrvolume = vector_dot(temp, v);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume == 0.0)
        return 0.0;
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    minsine2 = HUGEQUAL;     /* start with absurdly big value for sine */
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
                sine2 = 0.0;
            }
            
            /* update minimum sine */
            if (sine2 < minsine2)
            {
                minsine2 = sine2;
            }
        }
    }
    
    return sqrt(minsine2) * pyrvolume;
}


__device__ float minsine(float p1[3], float p2[3], float p3[3], float p4[3])
{
	float t[3], u[3], v[3]; /* tet vectors */
	float temp[3];
    float edgelength[3][4]; /* the lengths of each of the edges of the tet */
    float facenormal[4][3]; /* the normals of each face of the tet */
    float dx, dy, dz;       /* intermediate values of edge lengths */
    float facearea2[4];     /* areas of the faces of the tet */
    float pyrvolume;        /* volume of tetrahedron */
    float sine2, minsine2;  /* the sine (squared) of the dihedral angle */
    int i, j, k, l;          /* loop indices */
	float point[4][3];

	for (i = 0; i < 3; i++)
	{
		point[0][i] = p1[i];
		point[1][i] = p2[i];
		point[2][i] = p3[i];
		point[3][i] = p4[i];
	}
    
    /* calculate the volume*6 of the tetrahedron */
	vector_minus(point[1], point[0], t);
	vector_minus(point[2], point[0], u);
	vector_minus(point[3], point[0], v);
	vector_cross(t, u, temp);
    pyrvolume = vector_dot(temp, v);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume == 0.0)
        return 0.0;
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    minsine2 = HUGEQUAL;     /* start with absurdly big value for sine */
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
                sine2 = 0.0;
            }
            
            /* update minimum sine */
            if (sine2 < minsine2)
            {
                minsine2 = sine2;
            }
        }
    }
    
    return sqrt(minsine2) * pyrvolume;
}


__device__ float tetquality(float *points, int *meshtets, int tetcnt, int tetidx, int qualmeasure)
{
	if (tetidx < 0 || tetidx > tetcnt-1)
		return -1.0;

	float point[4][3];
	float quality = 0.0; /* the quality of this tetrahedron */
	int pidx;

	for (int j = 0; j < 4; j++)
	{
		pidx = meshtets[4*tetidx+j];
		point[j][0] = points[3*pidx];
		point[j][1] = points[3*pidx+1];
		point[j][2] = points[3*pidx+2];
	}

	switch (qualmeasure)
	{
	case CUDA_QUAL_MINSINE:
		quality = minsine(point);
		break;
	case CUDA_QUAL_BIASEDMINSINE:
		//quality = biasedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MEANSINE:
		//quality = meansine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINSINEANDEDGERATIO:
		//quality = minsineandedgeratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_RADIUSRATIO:
		//quality = radiusratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_VLRMS3RATIO:
		//quality = vlrms3ratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_WARPEDMINSINE:
		//quality = warpedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, false);
		break;
	case CUDA_QUAL_MAXANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, true);
		break;
	}
	return quality;
}

__device__ float tetquality(float p1[3], float p2[3], float p3[3], float p4[3], int qualmeasure)
{
	float quality = 0.0; /* the quality of this tetrahedron */
	switch (qualmeasure)
	{
	case CUDA_QUAL_MINSINE:
		quality = minsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_BIASEDMINSINE:
		//quality = biasedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MEANSINE:
		//quality = meansine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINSINEANDEDGERATIO:
		//quality = minsineandedgeratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_RADIUSRATIO:
		//quality = radiusratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_VLRMS3RATIO:
		//quality = vlrms3ratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_WARPEDMINSINE:
		//quality = warpedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, false);
		break;
	case CUDA_QUAL_MAXANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, true);
		break;
	}
	return quality;
}


__device__ float tetquality(float point[4][3], int qualmeasure)
{
	float quality = 0.0; /* the quality of this tetrahedron */
	switch (qualmeasure)
	{
	case CUDA_QUAL_MINSINE:
		quality = minsine(point);
		break;
	case CUDA_QUAL_BIASEDMINSINE:
		//quality = biasedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MEANSINE:
		//quality = meansine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINSINEANDEDGERATIO:
		//quality = minsineandedgeratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_RADIUSRATIO:
		//quality = radiusratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_VLRMS3RATIO:
		//quality = vlrms3ratio(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_WARPEDMINSINE:
		//quality = warpedminsine(p1, p2, p3, p4);
		break;
	case CUDA_QUAL_MINANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, false);
		break;
	case CUDA_QUAL_MAXANGLE:
		//quality = minmaxangle(p1, p2, p3, p4, true);
		break;
	}
	return quality;
}

__device__ float mintetquality(float *points, int *meshtets, int tetcnt, int qualmeasure)
{
	float qual, minqual;
	float tetpoint[4][3];
	int pidx;
	minqual = HUGEQUAL;
	for (int i = 0; i < tetcnt; i++)
	{
		// get tetra points
		for (int j = 0; j < 4; j++)
		{
			pidx = meshtets[4*i+j];
			tetpoint[j][0] = points[3*pidx];
			tetpoint[j][1] = points[3*pidx+1];
			tetpoint[j][2] = points[3*pidx+2];
		}
		// calculate tetra quality
		qual = tetquality(tetpoint, qualmeasure);
		// fetch the minimum quality value
		if (minqual > qual)
			minqual = qual;
	}
	return minqual;
}

__device__ float minstackquality(float *points, int *meshtets, int *tethandlestack, int stackcnt, int qualmeasure)
{
	float qual, minqual;
	float tetpoint[4][3];
	int tetidx;
	int pidx;

	minqual = HUGEQUAL;
	for (int i = 0; i < stackcnt; i++)
	{
		tetidx = tethandlestack[i];
		for (int j = 0; j < 4; j ++)
		{
			pidx = meshtets[4*tetidx+j];
			tetpoint[j][0] = points[3*pidx];
			tetpoint[j][1] = points[3*pidx+1];
			tetpoint[j][2] = points[3*pidx+2];
		}
		qual = tetquality(tetpoint, qualmeasure);
		if (minqual > qual)
			minqual = qual;
	}
	return minqual;
}

__device__ float minstackquality(float *tetrapoints, int tetracnt, int qualmeasure)
{
	float qual, minqual;
	float tetpoint[4][3];

	minqual = HUGEQUAL;
	for (int i = 0; i < tetracnt; i++)
	{
		for (int j = 0; j < 4; j ++)
		{
			tetpoint[j][0] = tetrapoints[3*(4*i+j)];
			tetpoint[j][1] = tetrapoints[3*(4*i+j)+1];
			tetpoint[j][2] = tetrapoints[3*(4*i+j)+2];
		}
		qual = tetquality(tetpoint, qualmeasure);
		if (minqual > qual)
			minqual = qual;
	}
	return minqual;
}

__global__ void cuda_mintetquality(float *points, int *meshtets, int tetcnt, int qualmeasure, float *tetqual)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	float tetpoint[4][3];
	float qual;
	int pidx;

	while(tid < tetcnt)
	{
		// if the tet is invalid
		if (meshtets[4*tid] == -1)
		{
			tetqual[tid] = 1.0;
			tid += offset;
			continue;
		}

		// get tetra points
		for (int j = 0; j < 4; j++)
		{
			pidx = meshtets[4*tid+j];
			tetpoint[j][0] = points[3*pidx];
			tetpoint[j][1] = points[3*pidx+1];
			tetpoint[j][2] = points[3*pidx+2];
		}
		// calculate tetra quality
		qual = tetquality(tetpoint, qualmeasure);
		tetqual[tid] = qual;

		tid += offset;
	}
}

void cutetramesh_quality(float *points, int pointcnt, cu_tetra *ctetra, int tetcnt, int qualmeasure, float &minqual)
{
	int *meshtets;
	meshtets = new int[4*tetcnt];
	for (int i = 0; i < tetcnt; i++)
	{
		meshtets[4*i]   = ctetra[i].v[0];
		meshtets[4*i+1] = ctetra[i].v[1];
		meshtets[4*i+2] = ctetra[i].v[2];
		meshtets[4*i+3] = ctetra[i].v[3];
	}

	float *dev_points;
	int *dev_meshtets;
	float *dev_tetqual;

	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
	cudaMalloc((void**)&dev_tetqual, tetcnt*sizeof(float));
	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);

	int blocks = imin(tetcnt, BlockPerGrid);
	cuda_mintetquality<<<blocks,(tetcnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, tetcnt, qualmeasure, dev_tetqual);

	float *tetqual;
	tetqual = new float[tetcnt];
	cudaMemcpy(tetqual, dev_tetqual, tetcnt*sizeof(float), cudaMemcpyDeviceToHost);

	minqual = HUGEQUAL;
	for (int i = 0; i < tetcnt; i++)
		if (tetqual[i] < minqual)
			minqual = tetqual[i];

	cudaFree(dev_points);
	cudaFree(dev_meshtets);
	cudaFree(dev_tetqual);
	delete [] tetqual;
	delete [] meshtets;
}

extern "C" void cuda_tetquality(float *points, int pointcnt, int *meshtets, int tetcnt, int qualmeasure, float &minqual)
{
	float *dev_points;
	int *dev_meshtets;
	float *dev_tetqual;

	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
	cudaMalloc((void**)&dev_tetqual, tetcnt*sizeof(float));
	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);

	int blocks = imin(tetcnt, BlockPerGrid);
	cuda_mintetquality<<<blocks,(tetcnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, tetcnt, qualmeasure, dev_tetqual);

	float *tetqual;
	tetqual = new float[tetcnt];
	cudaMemcpy(tetqual, dev_tetqual, tetcnt*sizeof(float), cudaMemcpyDeviceToHost);

	minqual = HUGEQUAL;
	for (int i = 0; i < tetcnt; i++)
		if (tetqual[i] > MINFLIPIMPROVE  && tetqual[i] < minqual)
			minqual = tetqual[i];

	cudaFree(dev_points);
	cudaFree(dev_meshtets);
	cudaFree(dev_tetqual);
	delete [] tetqual;
}
/********************* end of tetrahedron quality ************************/

/************************ Parallel vertex smoothing ***************************/
__global__ void laplacianSmoothing(float *points, int pointcnt, int *neighbour, int *neighbourcnt, int largestn,
								   int *incidenttet, int *incidenttetcnt, int largesttet,
								   int *meshtets, int tetcnt)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	float point[3];
	float oldpoint[3];
	float qualbefore, qualafter;

	while(tid < pointcnt)
	{
		// calculate original quality
		qualbefore = minstackquality(points, meshtets, &incidenttet[largesttet*tid], incidenttetcnt[tid], CUDA_QUAL_MINSINE);

		// fetch neighbour count
		int ncnt = neighbourcnt[tid];
		int nbr;
		point[0] = point[1] = point[2] = 0.0;

		// save original point coordinates
		oldpoint[0] = points[3*tid];
		oldpoint[1] = points[3*tid+1];
		oldpoint[2] = points[3*tid+2];

		// calculate new point coordinates
		for (int i = 0; i < ncnt; i++)
		{
			nbr = neighbour[tid*largestn+i];
			point[0] += points[3*nbr];
			point[1] += points[3*nbr+1];
			point[2] += points[3*nbr+2];
		}
		point[0] /= ncnt;
		point[1] /= ncnt;
		point[2] /= ncnt;

		// set new point
		points[3*tid]   = point[0];
		points[3*tid+1] = point[1];
		points[3*tid+2] = point[2];

		// calculate new tetra quality
		qualafter = minstackquality(points, meshtets, &incidenttet[largesttet*tid], incidenttetcnt[tid], CUDA_QUAL_MINSINE);
		if (qualafter < qualbefore)
		{
			points[3*tid]   = oldpoint[0];
			points[3*tid+1] = oldpoint[1];
			points[3*tid+2] = oldpoint[2];
		}

		tid += offset;
	}
}


// grouped smoothing：only smooth one group at one time
__global__ void group_laplacianSmoothing(float *points, int pointcnt, int *neighbour, int *neighbourcnt, int largestn,
										  int *pointgroup, int activegroup,int *incidenttet, int *incidenttetcnt, 
										  int largesttet, int *meshtets, int tetcnt)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	float point[3];
	float oldpoint[3];
	float qualbefore, qualafter;

	while(tid < pointcnt && pointgroup[tid] == activegroup)
	{
		// calculate original quality
		qualbefore = minstackquality(points, meshtets, &incidenttet[largesttet*tid], incidenttetcnt[tid], CUDA_QUAL_MINSINE);

		// fetch neighbour count
		int ncnt = neighbourcnt[tid];
		int nbr;
		point[0] = point[1] = point[2] = 0.0;

		// save original point coordinates
		oldpoint[0] = points[3*tid];
		oldpoint[1] = points[3*tid+1];
		oldpoint[2] = points[3*tid+2];

		// calculate new point coordinates
		for (int i = 0; i < ncnt; i++)
		{
			nbr = neighbour[tid*largestn+i];
			point[0] += points[3*nbr];
			point[1] += points[3*nbr+1];
			point[2] += points[3*nbr+2];
		}
		point[0] /= ncnt;
		point[1] /= ncnt;
		point[2] /= ncnt;

		// set new point
		points[3*tid]   = point[0];
		points[3*tid+1] = point[1];
		points[3*tid+2] = point[2];

		// calculate new tetra quality
		qualafter = minstackquality(points, meshtets, &incidenttet[largesttet*tid], incidenttetcnt[tid], CUDA_QUAL_MINSINE);
		if (qualafter < qualbefore)
		{
			points[3*tid]   = oldpoint[0];
			points[3*tid+1] = oldpoint[1];
			points[3*tid+2] = oldpoint[2];
		}

		tid += offset;
	}
}

extern "C" void cuda_vertexSmoothing(float *points, int pointcnt, int *hneighbour, int *neighbourcnt, int *pointgroup, int pointgroupcnt, 
									 float *newpoints, int *hincidenttet, int *incidenttetcnt, int *meshtets, int tetcnt, 
									 int largestn, int largesttet, int smoothpasscnt, float &time)
{
	int hncnt, hntet;
	hncnt = pointcnt * largestn;
	hntet = pointcnt * largesttet;

	// 分配设备存储空间
	int *dev_neighbour;
	int *dev_neighbourcnt;
	float *dev_points;
	int *dev_incidenttet;
	int *dev_incidenttetcnt;
	int *dev_meshtets;
	int *dev_pointgroup;
	bool *dev_pointsmoothed;

	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_neighbour, hncnt*sizeof(int));
	cudaMalloc((void**)&dev_neighbourcnt, pointcnt*sizeof(int));
	cudaMalloc((void**)&dev_incidenttet, hntet*sizeof(int));
	cudaMalloc((void**)&dev_incidenttetcnt, pointcnt*sizeof(int));
	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
	cudaMalloc((void**)&dev_pointgroup, pointcnt*sizeof(int));
	cudaMalloc((void**)&dev_pointsmoothed, pointcnt*sizeof(bool));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neighbour, hneighbour, pointcnt*largestn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neighbourcnt, neighbourcnt, pointcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_incidenttet, hincidenttet, hntet*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_incidenttetcnt, incidenttetcnt, pointcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pointgroup, pointgroup, pointcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemset(dev_pointsmoothed, 0, pointcnt*sizeof(bool));

	int blocks = imin(pointcnt, BlockPerGrid);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// 调用核函数
	//int blockPerGrid = (threadsPerBlock+pointcnt)/threadsPerBlock;
	while (smoothpasscnt)
	{
		for (int passidx = 0; passidx < smoothpasscnt; passidx ++)
		{
			//分 groupcnt 次发射 kernel
			for (int i = 0; i < pointgroupcnt; i++)
			{
				group_laplacianSmoothing<<<blocks, (pointcnt+blocks-1)/blocks>>>(dev_points, pointcnt, dev_neighbour, dev_neighbourcnt, largestn, dev_pointgroup, i, 
					dev_incidenttet, dev_incidenttetcnt, largesttet, dev_meshtets, tetcnt);
			}
		}
		--smoothpasscnt;
	}
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);
	time = elaspsedTime;

	cudaMemcpy(newpoints, dev_points, 3*pointcnt*sizeof(float), cudaMemcpyDeviceToHost);

	bool *pointsmoothed;
	pointsmoothed = new bool[pointcnt];
	cudaMemcpy(pointsmoothed, dev_pointsmoothed, pointcnt*sizeof(bool), cudaMemcpyDeviceToHost);


	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	cudaFree(dev_points);
	cudaFree(dev_neighbour);
	cudaFree(dev_neighbourcnt);
	cudaFree(dev_incidenttet);
	cudaFree(dev_incidenttetcnt);
	cudaFree(dev_meshtets);
	cudaFree(dev_pointgroup);
}

/************************** end of smoothing *******************************/

/************************** Parallel Flipping ******************************/
__global__ void flip23_explore(float *points, int *meshtets, struct cu_flip23face *face, int facecnt, 
							   struct cu_halfface *halfface, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int facepoint[3];
	int toppoint[2];
	int tet[2];
	int tetpoint[8];
	int newtet[12];
	float qualbefore, qualafter;
	//struct cu_flip23face currface;

	while(tid < facecnt)
	{
		//currface = face[tid];

		// check if it is a boundary face
		if (face[tid].hf[0] == -1 || face[tid].hf[1] == -1)
		{
			tid += offset;
			continue;
		}

		// get face points
		facepoint[0] = halfface[face[tid].hf[0]].pointhandle[0];
		facepoint[1] = halfface[face[tid].hf[0]].pointhandle[1];
		facepoint[2] = halfface[face[tid].hf[0]].pointhandle[2];

		// get tet relative data : tet[2]  tetpoint[2][4]  toppoint[2]
		for (int i = 0; i < 2; i ++)
		{
			// get two tets incident to the current face
			tet[i] = face[tid].hf[i]>>2;
			tetpoint[i<<2]     = meshtets[tet[i]<<2];
			tetpoint[(i<<2)+1] = meshtets[(tet[i]<<2)+1];
			tetpoint[(i<<2)+2] = meshtets[(tet[i]<<2)+2];
			tetpoint[(i<<2)+3] = meshtets[(tet[i]<<2)+3];

			for (int j = 0; j < 4; j++)
			{
				if (tetpoint[i*4+j] != facepoint[0] &&
					tetpoint[i*4+j] != facepoint[1] &&
					tetpoint[i*4+j] != facepoint[2])
				{
					toppoint[i] = tetpoint[i*4+j];
					break;
				}
			}
		}

		// get new tets
		newtet[0]  = toppoint[0]; 
		newtet[1]  = facepoint[0];
		newtet[2]  = toppoint[1];
		newtet[3]  = facepoint[2];

		newtet[4]  = toppoint[0];
		newtet[5]  = facepoint[0];
		newtet[6]  = facepoint[1];
		newtet[7]  = toppoint[1];

		newtet[8]  = toppoint[0];
		newtet[9]  = facepoint[1];
		newtet[10] = facepoint[2];
		newtet[11] = toppoint[1];

		// calculate original and new quality
		qualbefore = mintetquality(points, tetpoint, 2, qualmeasure);
		qualafter = mintetquality(points, newtet, 3, qualmeasure);

		// set face value
		face[tid].quality = qualbefore;
		face[tid].val = qualafter - qualbefore;

		tid += offset;
	}
}


__global__ void flip23 (int *meshtets, int tetcnt, struct cu_flip23face *face, int facecnt, struct cu_halfface *halfface, int *selectface, int selectfacecnt)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int facepoint[3];
	int toppoint[2];
	int tet[2];
	int tetpoint[8];
	int pidx;
	int tetpidx[2][5];   // record the tetra inner index of five points
	int newtet;
	cu_flip23face currface;
	int faceset[9];
	int i, j, k;

	while(tid < selectfacecnt)
	{
		currface = face[selectface[tid]];
		faceset[0] = selectface[tid];

		// get face points
		facepoint[0] = halfface[currface.hf[0]].pointhandle[0];
		facepoint[1] = halfface[currface.hf[0]].pointhandle[1];
		facepoint[2] = halfface[currface.hf[0]].pointhandle[2];

		// get tet relative data : tet[2]  tetpoint[2][4]  toppoint[2]
		for (i = 0; i < 2; i ++)
		{
			// get two tets incident to the current face
			tet[i] = currface.hf[i]>>2;
			tetpoint[i<<2]     = meshtets[tet[i]<<2];
			tetpoint[(i<<2)+1] = meshtets[(tet[i]<<2)+1];
			tetpoint[(i<<2)+2] = meshtets[(tet[i]<<2)+2];
			tetpoint[(i<<2)+3] = meshtets[(tet[i]<<2)+3];

			for (j = 0; j < 4; j++)
			{
				for (k = 0; k < 3; k++)
				{
					if (tetpoint[(i<<2)+j] == facepoint[k])
					{
						tetpidx[i][2+k] = j;
						break;
					}
				}

				if (tetpoint[i*4+j] != facepoint[0] &&
					tetpoint[i*4+j] != facepoint[1] &&
					tetpoint[i*4+j] != facepoint[2])
				{
					toppoint[i] = tetpoint[i*4+j];
					tetpidx[i][i] = j;
				}
			}
		}

		faceset[1] = halfface[(tet[0]<<2)+tetpidx[0][2]].face;
		faceset[2] = halfface[(tet[0]<<2)+tetpidx[0][3]].face;
		faceset[3] = halfface[(tet[0]<<2)+tetpidx[0][4]].face;
		faceset[4] = halfface[(tet[1]<<2)+tetpidx[1][2]].face;
		faceset[5] = halfface[(tet[1]<<2)+tetpidx[1][3]].face;
		faceset[6] = halfface[(tet[1]<<2)+tetpidx[1][4]].face;

		/* Set new tetras */
		// new tet0
		pidx = tet[0]<<2;
		meshtets[pidx]   = toppoint[0];   meshtets[pidx+1] = facepoint[0];
		meshtets[pidx+2] = toppoint[1];   meshtets[pidx+3] = facepoint[2];

		halfface[pidx].pointhandle[0] = meshtets[pidx+1];
		halfface[pidx].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+1].pointhandle[0] = meshtets[pidx];
		halfface[pidx+1].pointhandle[1] = meshtets[pidx+3];
		halfface[pidx+1].pointhandle[2] = meshtets[pidx+2];

		halfface[pidx+2].pointhandle[0] = meshtets[pidx];
		halfface[pidx+2].pointhandle[1] = meshtets[pidx+1];
		halfface[pidx+2].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+3].pointhandle[0] = meshtets[pidx];
		halfface[pidx+3].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx+3].pointhandle[2] = meshtets[pidx+1];

		// new tet1
		pidx = tet[1]<<2;
		meshtets[pidx]   = toppoint[0];   meshtets[pidx+1] = facepoint[0];
		meshtets[pidx+2] = facepoint[1];  meshtets[pidx+3] = toppoint[1];

		halfface[pidx].pointhandle[0] = meshtets[pidx+1];
		halfface[pidx].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+1].pointhandle[0] = meshtets[pidx];
		halfface[pidx+1].pointhandle[1] = meshtets[pidx+3];
		halfface[pidx+1].pointhandle[2] = meshtets[pidx+2];

		halfface[pidx+2].pointhandle[0] = meshtets[pidx];
		halfface[pidx+2].pointhandle[1] = meshtets[pidx+1];
		halfface[pidx+2].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+3].pointhandle[0] = meshtets[pidx];
		halfface[pidx+3].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx+3].pointhandle[2] = meshtets[pidx+1];

		// add new tet
		newtet = tid+tetcnt;
		pidx = newtet<<2;
		meshtets[pidx]   = toppoint[0];  meshtets[pidx+1] = facepoint[1];
		meshtets[pidx+2] = facepoint[2]; meshtets[pidx+3] = toppoint[1];

		halfface[pidx].pointhandle[0] = meshtets[pidx+1];
		halfface[pidx].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+1].pointhandle[0] = meshtets[pidx];
		halfface[pidx+1].pointhandle[1] = meshtets[pidx+3];
		halfface[pidx+1].pointhandle[2] = meshtets[pidx+2];

		halfface[pidx+2].pointhandle[0] = meshtets[pidx];
		halfface[pidx+2].pointhandle[1] = meshtets[pidx+1];
		halfface[pidx+2].pointhandle[2] = meshtets[pidx+3];

		halfface[pidx+3].pointhandle[0] = meshtets[pidx];
		halfface[pidx+3].pointhandle[1] = meshtets[pidx+2];
		halfface[pidx+3].pointhandle[2] = meshtets[pidx+1];

		// add two new face
		face[(facecnt+tid)*2].hf[0] = 4*tet[0]+1;
		face[(facecnt+tid)*2].hf[1] = 4*newtet+1;
		face[(facecnt+tid)*2+1].hf[0] = 4*tet[1]+1;
		face[(facecnt+tid)*2+1].hf[1] = 4*newtet+2;

		faceset[7] = (facecnt+tid)*2;
		faceset[8] = (facecnt+tid)*2+1;

		halfface[4*tet[0]+1].face = (facecnt+tid)*2;
		halfface[4*newtet+1].face = (facecnt+tid)*2;
		halfface[4*tet[1]+1].face = (facecnt+tid)*2+1;
		halfface[4*newtet+2].face = (facecnt+tid)*2+1;

		// update some old face and halfface relationship
		face[selectface[tid]].hf[0] = 4*tet[1]+2;
		face[selectface[tid]].hf[1] = 4*tet[0]+3;
		halfface[4*tet[1]+2].face = selectface[tid];
		halfface[4*tet[0]+3].face = selectface[tid];

		// update new tet0
		halfface[4*tet[0]].face = faceset[5];
		if (face[faceset[5]].hf[0]>>2 == tet[1])
			face[faceset[5]].hf[0] = 4*tet[0];
		else
			face[faceset[5]].hf[1] = 4*tet[0];

		halfface[4*tet[0]+2].face = faceset[2];
		if (face[faceset[2]].hf[0]>>2 == tet[0])
			face[faceset[2]].hf[0] = 4*tet[0]+2;
		else
			face[faceset[2]].hf[1] = 4*tet[0]+2;

		// update new tet1
		halfface[4*tet[1]].face = faceset[6];
		if (face[faceset[6]].hf[0]>>2 == tet[1])
			face[faceset[6]].hf[0] = 4*tet[1];
		else
			face[faceset[6]].hf[1] = 4*tet[1];

		halfface[4*tet[1]+3].face = faceset[3];
		if (face[faceset[3]].hf[0]>>2 == tet[0])
			face[faceset[3]].hf[0] = 4*tet[1]+3;
		else
			face[faceset[3]].hf[1] = 4*tet[1]+3;

		// update newtet
		halfface[4*newtet].face = faceset[4];
		if (face[faceset[4]].hf[0]>>2 == tet[1])
			face[faceset[4]].hf[0] = 4*newtet;
		else
			face[faceset[4]].hf[1] = 4*newtet;

		halfface[4*newtet+3].face = faceset[1];
		if (face[faceset[1]].hf[0]>>2 == tet[0])
			face[faceset[1]].hf[0] = 4*newtet + 3;
		else
			face[faceset[1]].hf[1] = 4*newtet + 3;

		// add offset and go next loop
		tid += offset;
	}
}

void faceSelecting(struct cu_flip23face *face, int facecnt, int tetcnt, int* &newface, int &newfacecnt)
{
	int selectfacecnt = 0;
	int* selectface = new int[facecnt];
	int i, j;
	float threashold = float(1.0e-5);

	/* pick out the faces by flipping succeed val
	and sort it*/
	for (i = 0; i < facecnt; i++)
	{
		if (!(face[i].val > threashold))
			continue;

		// insert the face
		for (j = selectfacecnt-1; j > -1 ; j--)
		{
			if (face[selectface[j]].quality > face[i].quality)
				selectface[j+1] = selectface[j];
			else 
				break;
		}
		selectface[j+1] = i;
		++ selectfacecnt;
	}

	/* if there are no any face meeting the requirement */
	if (!selectfacecnt)
	{
		newfacecnt = 0;
		newface = NULL;
		return;
	}

	/* select a new face set in which any two of them not in a same tet */
	int tet[2];
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, 1, tetcnt*sizeof(bool));

	// push the first face
	newface = new int[facecnt];
	newface[0] = selectface[0];
	newfacecnt = 1;

	// set flags of tets incident to the face
	tetflag[(face[newface[0]].hf[0])>>2] = 0;
	tetflag[(face[newface[0]].hf[1])>>2] = 0;

	/* if the tets incident to a face are available, 
	then add the face into the array*/
	for (i = 1; i < selectfacecnt; i++)
	{
		tet[0] = (face[selectface[i]].hf[0])>>2;
		tet[1] = (face[selectface[i]].hf[1])>>2;
		if (tetflag[tet[0]] && tetflag[tet[1]])
		{
			newface[newfacecnt++] = selectface[i];
			tetflag[tet[0]] = 0;
			tetflag[tet[1]] = 0;
		}
	}
}

extern "C" void cuda_flip23(float *points, int pointcnt, int *meshtets_, int tetcnt, int *face, int facecnt, 
							int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
							int &flipsucc, float &time)
{
	// create face and halfface structure
	int *meshtets, *tempmeshtets;
	struct cu_flip23face *cface, *tempcafce;
	struct cu_halfface *chalfface, *tempchalfface;
	float qualbefore;
	float qualafter;
	int tetcapacity;
	int facecapacity;
	int halffacecapacity;

	tetcapacity = int(1.1*tetcnt);
	facecapacity = int(1.1*facecnt);
	halffacecapacity = int(1.1*halffacecnt);
	meshtets = new int[4*tetcapacity];
	cface = new struct cu_flip23face[facecapacity];
	chalfface = new struct cu_halfface[halffacecapacity];

	memcpy(meshtets, meshtets_, 4*tetcnt*sizeof(int));

	// calculate quality before flipping
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

	// CUDA Parallel
	// 分配设备存储空间
	int loop = 1;
	float *dev_points;
	int *dev_meshtets, *dev_tempmeshtets;
	struct cu_flip23face *dev_face, *dev_tempface;
	struct cu_halfface *dev_halfface, *dev_temphalfface;
	int *selectface;
	int selectfacecnt;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (int i = 0; i < facecnt; i++)
	{
		cface[i].hf[0] = face[2*i];
		cface[i].hf[1] = face[2*i+1];
		cface[i].quality = -1;
		cface[i].val = -1;
	}

	for (int i = 0; i < halffacecnt; i++)
	{
		chalfface[i].pointhandle[0] = halfface[4*i];
		chalfface[i].pointhandle[1] = halfface[4*i+1];
		chalfface[i].pointhandle[2] = halfface[4*i+2];
		chalfface[i].face = halfface[4*i+3];
	}


	// flip23 会生成新的四面体、半面和面，因此先预分配一部分空间
	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_meshtets, 4*tetcapacity*sizeof(int));
	cudaMalloc((void**)&dev_face, facecapacity*sizeof(cu_flip23face));
	cudaMalloc((void**)&dev_halfface, halffacecapacity*sizeof(cu_halfface));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_face, cface, facecnt*sizeof(cu_flip23face), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_halfface, chalfface, halffacecnt*sizeof(cu_halfface), cudaMemcpyHostToDevice);

	while(loop)
	{
		// 调用核函数
		int blocks = imin(pointcnt, BlockPerGrid);

		/* strategy : flip test -> select proper face set -> do flipping */
		// step1 : flip test
		flip23_explore<<<blocks, (facecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_face, facecnt, dev_halfface, qualmeasure);
		cudaMemcpy(cface, dev_face, facecnt*sizeof(cu_flip23face), cudaMemcpyDeviceToHost);

		// step2 : select proper face set
		faceSelecting(cface, facecnt, tetcnt, selectface, selectfacecnt);

		// step3 : do flipping
		if (selectfacecnt)
		{
			// check if the memory is not enough, apply more
			if (tetcnt + selectfacecnt > tetcapacity)
			{
				tetcapacity = int(1.2*tetcapacity);
				cudaMalloc((void**)&dev_tempmeshtets, 4*tetcnt*sizeof(int));
				cudaMemcpy(dev_tempmeshtets, dev_meshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToDevice);
				cudaFree(dev_meshtets);
				cudaMalloc((void**)&dev_meshtets, 4*tetcapacity*sizeof(int));
				cudaMemcpy(dev_meshtets, dev_tempmeshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToDevice);
				cudaFree(dev_tempmeshtets);

				tempmeshtets = new int[4*tetcnt];
				memcpy(tempmeshtets, meshtets, 4*tetcnt*sizeof(int));
				delete [] meshtets;
				meshtets = new int[4*tetcapacity];
				memcpy(meshtets, tempmeshtets, 4*tetcnt*sizeof(int));
				delete [] tempmeshtets;

				halffacecapacity = int(1.2*halffacecapacity);
				cudaMalloc((void**)&dev_temphalfface, halffacecnt*sizeof(struct cu_halfface));
				cudaMemcpy(dev_temphalfface, dev_halfface, halffacecnt*sizeof(struct cu_halfface), cudaMemcpyDeviceToDevice);
				cudaFree(dev_halfface);
				cudaMalloc((void**)&dev_halfface, halffacecapacity*sizeof(struct cu_halfface));
				cudaMemcpy(dev_halfface, dev_temphalfface, halffacecnt*sizeof(struct cu_halfface), cudaMemcpyDeviceToDevice);
				cudaFree(dev_temphalfface);

				tempchalfface = new struct cu_halfface[halffacecnt];
				memcpy(tempchalfface, chalfface, halffacecnt*sizeof(cu_halfface));
				delete [] chalfface;
				chalfface = new struct cu_halfface[halffacecapacity];
				memcpy(chalfface, tempchalfface, halffacecnt*sizeof(cu_halfface));
				delete [] tempchalfface;
			}
			if (facecnt + 3*selectfacecnt > facecapacity)
			{
				facecapacity = int(1.2*facecapacity);
				cudaMalloc((void**)&dev_tempface, facecnt*sizeof(struct cu_flip23face));
				cudaMemcpy(dev_tempface, dev_face, facecnt*sizeof(struct cu_flip23face), cudaMemcpyDeviceToDevice);
				cudaFree(dev_face);
				cudaMalloc((void**)&dev_face, facecapacity*sizeof(struct cu_flip23face));
				cudaMemcpy(dev_face, dev_tempface, facecnt*sizeof(struct cu_flip23face), cudaMemcpyDeviceToDevice);
				cudaFree(dev_tempface);

				tempcafce = new struct cu_flip23face[facecnt];
				memcpy(tempcafce, cface, facecnt*sizeof(struct cu_flip23face));
				delete [] cface;
				cface = new struct cu_flip23face[facecapacity];
				memcpy(cface, tempcafce, facecnt*sizeof(struct cu_flip23face));
				delete [] tempcafce;
			}
 
			int *dev_selectface;
			cudaMalloc((void**)&dev_selectface, selectfacecnt*sizeof(int));
			cudaMemcpy(dev_selectface, selectface, selectfacecnt*sizeof(int), cudaMemcpyHostToDevice);
			flip23<<<blocks, (selectfacecnt+blocks-1)/blocks>>>(dev_meshtets, tetcnt, dev_face, facecnt, dev_halfface, dev_selectface, selectfacecnt);

			tetcnt += selectfacecnt;
			facecnt += 2*selectfacecnt;
			halffacecnt += 4*selectfacecnt;

			cudaMemcpy(meshtets, dev_meshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToHost);
			cudaFree(dev_selectface);
		}
		/* end of strategy2 */

		// malloc in faceSelecting function
		delete [] selectface;
		--loop;
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);

	time = elaspsedTime;

	// calculate quality after flipping
	qualafter = 1.0;
	cudaMemcpy(meshtets, dev_meshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToHost);
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualafter);

	qualbefore_ = qualbefore;
	qualafter_ = qualafter;
	flipsucc = selectfacecnt;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_points);
	cudaFree(dev_meshtets);
	cudaFree(dev_halfface);
	cudaFree(dev_face);

	delete [] cface;
	delete [] chalfface;
}

__global__ void flip32_explore(float *points, int *meshtets, struct cu_flip32edge *edge, int edgecnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int tet[12];
	int newtet[8];
	int pidxvec[5];
	//int idx;
	int i,j;
	//int tetidx;
	//int pidx;
	float qualbefore, qualafter;
	//cu_flip32edge curredge;

	while(tid < edgecnt)
	{
		//curredge = edge[tid];

		/* fetch local data */
		// edge incident tets
		for (i = 0; i < 3; i ++)
		{
			//tetidx = edge[tid].tet[i];
			for (j = 0; j < 4; j++)
				tet[i*4+j] = meshtets[edge[tid].tet[i]*4+j];
		}

		/* flip32 incident five points: 
		   two are the endpoints of edge, 
		   the other three have to get from incident tets*/
		pidxvec[0] = edge[tid].p[0];
		pidxvec[1] = edge[tid].p[1];

		j = 2;
		// get two other points from one tet
		for (i = 0; i < 4; i++)
		{
			//pidx = tet[i];
			if (tet[i] != pidxvec[0] && tet[i] != pidxvec[1])
			{
				pidxvec[j] = tet[i];
				++j;
			}
		}

		// wrong data
		if (j != 4)
		{
			tid += offset;
			continue;
		}

		// get the last one from another tet
		pidxvec[4] = -1;
		for (i = 4; i < 12; i++)
		{
			//pidx = tet[i+4];
			if (tet[i] != pidxvec[0] && tet[i] != pidxvec[1] &&
				tet[i] != pidxvec[2] && tet[i] != pidxvec[3])
			{
				pidxvec[4] = tet[i];
				break;
			}
		}

		if (pidxvec[4] == -1)
		{
			tid += offset;
			continue;
		}
		/* end of data fetch */

		// calculate quality before flip32
		qualbefore = mintetquality(points, tet, 3, qualmeasure);

		// get new tets : the five points are not in order, so we probably try twice
		// conbination 1
		newtet[0] = pidxvec[0];
		newtet[1] = pidxvec[2];
		newtet[2] = pidxvec[3];
		newtet[3] = pidxvec[4];

		newtet[4] = pidxvec[1];
		newtet[5] = pidxvec[2];
		newtet[6] = pidxvec[4];
		newtet[7] = pidxvec[3];

		edge[tid].order = 0;

		qualafter = mintetquality(points, newtet, 2, qualmeasure);

		// if generate a reverse tet, try another way
		if (qualafter < 0)
		{
			newtet[0] = pidxvec[0];
			newtet[1] = pidxvec[2];
			newtet[2] = pidxvec[4];
			newtet[3] = pidxvec[3];

			newtet[4] = pidxvec[1];
			newtet[5] = pidxvec[2];
			newtet[6] = pidxvec[3];
			newtet[7] = pidxvec[4];

			edge[tid].order = 1;

			qualafter = mintetquality(points, newtet, 2, qualmeasure);
		}

		// set data
		edge[tid].quality = qualbefore;
		edge[tid].val = qualafter - qualbefore;

		tid += offset;
	}
}

__global__ void flip32(float *points, int *meshtets, struct cu_flip32edge *edge, 
					   int *selectedge, int selectedgecnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int tet[12];
	int newtet[8];
	int pidxvec[5];
	int i, j;
	//int pidx;
	//float qual;
	cu_flip32edge curredge;

	while(tid < selectedgecnt)
	{
		curredge = edge[selectedge[tid]];

		/* fetch local data */
		// edge incident tets
		for (i = 0; i < 3; i ++)
		{
			for (j = 0; j < 4; j++)
				tet[i*4+j] = meshtets[curredge.tet[i]*4+j];
		}

		/* flip32 incident five points: 
		   two are the endpoints of edge, 
		   the other three have to get from incident tets*/
		pidxvec[0] = curredge.p[0];
		pidxvec[1] = curredge.p[1];

		j = 2;
		// get two other points from one tet
		for (i = 0; i < 4; i++)
		{
			//pidx = tet[i];
			if (tet[i] != pidxvec[0] && tet[i] != pidxvec[1])
			{
				pidxvec[j] = tet[i];
				++j;
			}
		}

		// wrong data
		if (j != 4)
		{
			tid += offset;
			continue;
		}

		// get the last one from another tet
		pidxvec[4] = -1;
		for (i = 4; i < 12; i++)
		{
			//pidx = tet[i+4];
			if (tet[i] != pidxvec[0] && tet[i] != pidxvec[1] &&
				tet[i] != pidxvec[2] && tet[i] != pidxvec[3])
			{
				pidxvec[4] = tet[i];
				break;
			}
		}

		if (pidxvec[4] == -1)
		{
			tid += offset;
			continue;
		}
		/* end of data fetch */

		// calculate quality before flip32
		//qualbefore = mintetquality(points, tet, 3, qualmeasure);

		// get new tets : the five points are not in order, so we probably try twice
		// conbination 1
		if (curredge.order == 0)
		{
			newtet[0] = pidxvec[0];
			newtet[1] = pidxvec[2];
			newtet[2] = pidxvec[3];
			newtet[3] = pidxvec[4];

			newtet[4] = pidxvec[1];
			newtet[5] = pidxvec[2];
			newtet[6] = pidxvec[4];
			newtet[7] = pidxvec[3];
		}
		else
		{
			newtet[0] = pidxvec[0];
			newtet[1] = pidxvec[2];
			newtet[2] = pidxvec[4];
			newtet[3] = pidxvec[3];

			newtet[4] = pidxvec[1];
			newtet[5] = pidxvec[2];
			newtet[6] = pidxvec[3];
			newtet[7] = pidxvec[4];
		}


		// update data
		//if (qualafter - qualbefore > 1.0e-5)
		//{
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 4; j++)
				{
					meshtets[curredge.tet[i]*4+j] = newtet[i*4+j];
				}
			}

			for (i = 0; i < 4; i++)
			{
				meshtets[curredge.tet[2]*4+i] = -1;
			}
		//}
		tid += offset;
	}
}


void edgeSelecting(cu_flip32edge *cuedge, int edgecnt, int tetcnt, int *selectedge, int &selectedgecnt)
{
	int sedgecnt = 0;
	int* sedge = new int[edgecnt];
	int i, j;
	float threashold = float(1.0e-5);

	/* pick out the faces by flipping succeed val
	and sort it*/
	for (i = 0; i < edgecnt; i++)
	{
		if (!(cuedge[i].val > threashold))
			continue;

		// insert the face
		for (j = sedgecnt-1; j > -1 ; j--)
		{
			if (cuedge[sedge[j]].quality > cuedge[i].quality)
				sedge[j+1] = sedge[j];
			else 
				break;
		}
		sedge[j+1] = i;
		++ sedgecnt;
	}

	/* if there are no any face meeting the requirement */
	if (!sedgecnt)
	{
		selectedgecnt = 0;
		selectedge = NULL;
		return;
	}

	/* select a new face set in which any two of them not in a same tet */
	int tet[3];
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, 1, tetcnt*sizeof(bool));

	// push the first face
	selectedge[0] = sedge[0];
	selectedgecnt = 1;

	// set flags of tets incident to the face
	tetflag[cuedge[selectedge[0]].tet[0]] = 0;
	tetflag[cuedge[selectedge[0]].tet[1]] = 0;
	tetflag[cuedge[selectedge[0]].tet[2]] = 0;

	/* if the tets incident to a face are available, 
	then add the face into the array*/
	for (i = 1; i < sedgecnt; i++)
	{
		tet[0] = cuedge[sedge[i]].tet[0];
		tet[1] = cuedge[sedge[i]].tet[1];
		tet[2] = cuedge[sedge[i]].tet[2];
		if (tetflag[tet[0]] && tetflag[tet[1]] && tetflag[tet[2]])
		{
			selectedge[selectedgecnt++] = sedge[i];
			tetflag[tet[0]] = 0;
			tetflag[tet[1]] = 0;
			tetflag[tet[2]] = 0;
		}
	}
}

extern "C" void cuda_flip32(float *points, int pointcnt, int *meshtets, int tetcnt, int *flipedge, int edgecnt, 
				            int qualmeasure, float& qualbefore_, float& qualafter_, int &flipsucc, float &time)
{
	

	float qualbefore, qualafter;
	float *dev_points;
	int *dev_meshtets;
	struct cu_flip32edge *dev_edge;
	int *selectedge;
	int selectedgecnt;
	int loop = 4;

	// calculate quality before flip32
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	
	struct cu_flip32edge *cuedge;
	cuedge = new struct cu_flip32edge[edgecnt];

	for (int i = 0; i < edgecnt; i++)
	{
		cuedge[i].p[0] = flipedge[i*5];
		cuedge[i].p[1] = flipedge[i*5+1];
		cuedge[i].tet[0] = flipedge[i*5+2];
		cuedge[i].tet[1] = flipedge[i*5+3];
		cuedge[i].tet[2] = flipedge[i*5+4];
		cuedge[i].quality = 1.0;
		cuedge[i].val = -1.0;
	}

	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
	cudaMalloc((void**)&dev_edge, edgecnt*sizeof(cu_flip32edge));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_edge, cuedge, edgecnt*sizeof(cu_flip32edge), cudaMemcpyHostToDevice);
	
	flipsucc = 0;
	//while(loop)
	{
		/***********do parallel flip32**********/
		// flip32 explore
		int blocks = imin(edgecnt, BlockPerGrid);
		int threads = imin((edgecnt+blocks-1)/blocks, ThreadPerBlock);
		flip32_explore<<<blocks, threads>>>(dev_points, dev_meshtets, dev_edge, edgecnt, qualmeasure);
		cudaMemcpy(cuedge, dev_edge, edgecnt*sizeof(cu_flip32edge), cudaMemcpyDeviceToHost);

		// edge selecting
		selectedge = new int[edgecnt];
		edgeSelecting(cuedge, edgecnt, tetcnt, selectedge, selectedgecnt);

		// flip32 (set -1 to the meshtets[] of invalid tets)
		if (selectedgecnt)
		{
			int *dev_selectedge;
			cudaMalloc((void**)&dev_selectedge, selectedgecnt*sizeof(int));
			cudaMemcpy(dev_selectedge, selectedge, selectedgecnt*sizeof(int), cudaMemcpyHostToDevice);

			// do flip32
			blocks = imin(selectedgecnt, BlockPerGrid);
			flip32<<<blocks, (selectedgecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_edge, dev_selectedge, selectedgecnt, qualmeasure);

			flipsucc += selectedgecnt;
			cudaFree(dev_selectedge);
		}
		--loop;
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);

	time = elaspsedTime;

	// calculate quality after flipping
	cudaMemcpy(meshtets, dev_meshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToHost);
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualafter);

	qualbefore_ = qualbefore;
	qualafter_ = qualafter;

	delete [] selectedge;  // malloc in edge selecting function
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_points);
	cudaFree(dev_meshtets);
}

__global__ void flip32_explore_new(float *points, int *meshtets, cu_edge *edge, int edgecnt, cu_halfedge *halfedge, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int tet[12];
	int newtet[8];
	int pidxvec[5];
	int idx;
	int tetidx;
	int pidx;
	float qualbefore, qualafter;
	cu_edge curredge;

	while(tid < edgecnt)
	{
		curredge = edge[tid];

		// if the edge is boundary or the amount of incident tetras is not 3, go next
		if (curredge.is_boundary || curredge.halfedgecnt != 6)
		{
			tid += offset;
			continue;
		}

		/* fetch local data */
		// edge incident tets
		for (int i = 0; i < 3; i ++)
		{
			tetidx = curredge.halfedge[i*2]>>2;
			for (int j = 0; j < 4; j++)
				tet[i*4+j] = meshtets[tetidx*4+j];
		}

		/* flip32 incident five points: 
		   two are the endpoints of edge, 
		   the other three have to get from incident tets*/
		pidxvec[0] = halfedge[curredge.halfedge[0]].fromv;
		pidxvec[1] = halfedge[curredge.halfedge[0]].tov;

		idx = 0;
		// get two other points from one tet
		for (int i = 0; i < 4; i++)
		{
			pidx = tet[i];
			if (pidx != pidxvec[0] && pidx != pidxvec[1])
			{
				pidxvec[idx+2] = pidx;
				++idx;
			}
		}

		// wrong data
		if (idx != 2)
		{
			tid += offset;
			continue;
		}

		// get the last one from another tet
		for (int i = 0; i < 4; i++)
		{
			pidx = tet[i+4];
			if (pidx != pidxvec[0] && pidx != pidxvec[1] &&
				pidx != pidxvec[2] && pidx != pidxvec[3])
			{
				pidxvec[4] = pidx;
				break;
			}
		}
		/* end of data fetch */

		// calculate quality before flip32
		qualbefore = mintetquality(points, tet, 3, qualmeasure);

		// get new tets : the five points are not in order, so we probably try twice
		// conbination 1
		newtet[0] = pidxvec[0];
		newtet[1] = pidxvec[2];
		newtet[2] = pidxvec[3];
		newtet[3] = pidxvec[4];

		newtet[4] = pidxvec[1];
		newtet[5] = pidxvec[2];
		newtet[6] = pidxvec[4];
		newtet[7] = pidxvec[3];

		qualafter = mintetquality(points, newtet, 2, qualmeasure);

		// if generate a reverse tet, try another way
		if (qualafter < 0)
		{
			newtet[0] = pidxvec[0];
			newtet[1] = pidxvec[2];
			newtet[2] = pidxvec[4];
			newtet[3] = pidxvec[3];

			newtet[4] = pidxvec[1];
			newtet[5] = pidxvec[2];
			newtet[6] = pidxvec[3];
			newtet[7] = pidxvec[4];

			qualafter = mintetquality(points, newtet, 2, qualmeasure);
		}

		// set data
		edge[tid].quality = qualbefore;
		edge[tid].val = qualafter - qualbefore;

		tid += offset;
	}
}

void edgeSelecting_new(cu_edge *cuedge, int edgecnt, int tetcnt, int *selectedge, int &selectedgecnt)
{
	int sedgecnt = 0;
	int* sedge = new int[edgecnt];
	int i, j;
	float threashold = MINIMPROVEMENT;

	/* pick out the faces by flipping succeed val
	and sort it*/
	for (i = 0; i < edgecnt; i++)
	{
		if (!(cuedge[i].val > threashold))
			continue;

		// insert the face
		for (j = sedgecnt-1; j > -1 ; j--)
		{
			if (cuedge[sedge[j]].quality > cuedge[i].quality)
				sedge[j+1] = sedge[j];
			else 
				break;
		}
		sedge[j+1] = i;
		++ sedgecnt;
	}

	/* if there are no any face meeting the requirement */
	if (!sedgecnt)
	{
		selectedgecnt = 0;
		selectedge = NULL;
		return;
	}

	/* select a new face set in which any two of them not in a same tet */
	int tet[3];
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, 1, tetcnt*sizeof(bool));

	// push the first face
	selectedge[0] = sedge[0];
	selectedgecnt = 1;

	// set flags of tets incident to the face
	tetflag[cuedge[selectedge[0]].halfedge[0]>>2] = 0;
	tetflag[cuedge[selectedge[0]].halfedge[2]>>2] = 0;
	tetflag[cuedge[selectedge[0]].halfedge[4]>>2] = 0;

	/* if the tets incident to a face are available, 
	then add the face into the array*/
	for (i = 1; i < sedgecnt; i++)
	{
		tet[0] = cuedge[sedge[i]].halfedge[0]>>2;
		tet[1] = cuedge[sedge[i]].halfedge[2]>>2;
		tet[2] = cuedge[sedge[i]].halfedge[4]>>2;
		if (tetflag[tet[0]] && tetflag[tet[1]] && tetflag[tet[2]])
		{
			selectedge[selectedgecnt++] = sedge[i];
			tetflag[tet[0]] = 0;
			tetflag[tet[1]] = 0;
			tetflag[tet[2]] = 0;
		}
	}
}

__global__ void flip32_new(float *points, int *meshtets, cu_edge *edge, cu_halfedge *halfedge, 
					       int *selectedge, int selectedgecnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int tet[12];
	int newtet[8];
	int pidxvec[5];
	int idx;
	int tetidx;
	int pidx;
	float qualbefore, qualafter;
	cu_edge curredge;

	while(tid < selectedgecnt)
	{
		curredge = edge[selectedge[tid]];

		/* fetch local data */
		// edge incident tets
		for (int i = 0; i < 3; i ++)
		{
			tetidx = curredge.halfedge[i*2]>>2;
			for (int j = 0; j < 4; j++)
				tet[i*4+j] = meshtets[tetidx*4+j];
		}

		/* flip32 incident five points: 
		   two are the endpoints of edge, 
		   the other three have to get from incident tets*/
		pidxvec[0] = halfedge[curredge.halfedge[0]].fromv;
		pidxvec[1] = halfedge[curredge.halfedge[0]].tov;

		idx = 0;
		// get two other points from one tet
		for (int i = 0; i < 4; i++)
		{
			pidx = tet[i];
			if (pidx != pidxvec[0] && pidx != pidxvec[1])
			{
				pidxvec[idx+2] = pidx;
				++idx;
			}
		}

		// wrong data
		if (idx != 2)
		{
			tid += offset;
			continue;
		}

		// get the last one from another tet
		for (int i = 0; i < 4; i++)
		{
			pidx = tet[i+4];
			if (pidx != pidxvec[0] && pidx != pidxvec[1] &&
				pidx != pidxvec[2] && pidx != pidxvec[3])
			{
				pidxvec[4] = pidx;
				break;
			}
		}
		/* end of data fetch */

		// calculate quality before flip32
		qualbefore = mintetquality(points, tet, 3, qualmeasure);

		// get new tets : the five points are not in order, so we probably try twice
		// conbination 1
		newtet[0] = pidxvec[0];
		newtet[1] = pidxvec[2];
		newtet[2] = pidxvec[3];
		newtet[3] = pidxvec[4];

		newtet[4] = pidxvec[1];
		newtet[5] = pidxvec[2];
		newtet[6] = pidxvec[4];
		newtet[7] = pidxvec[3];

		qualafter = mintetquality(points, newtet, 2, qualmeasure);

		// if generate a reverse tet, try another way
		if (qualafter < 0)
		{
			newtet[0] = pidxvec[0];
			newtet[1] = pidxvec[2];
			newtet[2] = pidxvec[4];
			newtet[3] = pidxvec[3];

			newtet[4] = pidxvec[1];
			newtet[5] = pidxvec[2];
			newtet[6] = pidxvec[3];
			newtet[7] = pidxvec[4];

			qualafter = mintetquality(points, newtet, 2, qualmeasure);
		}

		// update data
		if (qualafter - qualbefore > 1.0e-5)
		{
			// update tetras
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					meshtets[(curredge.halfedge[i*2]>>2)<<2+j] = newtet[i*4+j];

					// update halfface
				}
			}

			for (int i = 0; i < 4; i++)
			{
				meshtets[(curredge.halfedge[4]>>2)<<2+i] = -1;
			}

			// update edges
			curredge.halfedgecnt = 0;  // current edge has been deleted

			//for (int i = 0; i < 3; i++)
			//{
			//	cuPrintf("Tetra %d : ", curredge.tet[i]);
			//	for (int j = 0; j < 4; j++)
			//	{
			//		cuPrintf("%d ", meshtets[curredge.tet[i]*4+j]);
			//	}
			//	cuPrintf("\n");
			//}
			//cuPrintf("\nqualbefore: %f   qualafter: %f \n", qualbefore, qualafter);
		}

		tid += offset;
	}
}

extern "C" void cuda_flip32_new(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
								int *halfedge, int halfedgecnt, int qualmeasure, float &qualbefore_, float &qualafeter_,
								int &flipsucc, float &time)
{
	cu_edge *cedge;
	cu_halfedge *chalfedge;
	cedge = new cu_edge[edgecnt];
	chalfedge = new cu_halfedge[halfedgecnt];

	// build edge
	int idx = 0;
	for (int i = 0; i < edgecnt; i ++)
	{
		cedge[i].is_boundary = edge[idx++];
		cedge[i].halfedgecnt = edge[idx++];
		for (int j = 0; j < cedge[i].halfedgecnt; j++)
		{
			cedge[i].halfedge[j] = edge[idx++];
		}
		cedge[i].quality = cedge[i].val = -1;
	}

	// build halfedge
	for (int i = 0; i < halfedgecnt; i++)
	{
		chalfedge[i].edgehandle = halfedge[i*3];
		chalfedge[i].fromv = halfedge[i*3+1];
		chalfedge[i].tov = halfedge[i*3+2];
	}

	float qualbefore;
	//float qualafter;
	float *dev_points;
	int *dev_meshtets;
	struct cu_edge *dev_edge;
	struct cu_halfedge *dev_halfedge;
	int *selectedge;
	int selectedgecnt;
	int loop = 1;

	// calculate quality before flip32
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
	cudaMalloc((void**)&dev_edge, edgecnt*sizeof(cu_edge));
	cudaMalloc((void**)&dev_halfedge, halfedgecnt*sizeof(cu_halfedge));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_edge, cedge, edgecnt*sizeof(cu_edge), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_halfedge, chalfedge, halfedgecnt*sizeof(cu_halfedge), cudaMemcpyHostToDevice);

	flipsucc = 0;
	while(loop)
	{

		/***********do parallel flip32**********/
		// flip32 explore
		int blocks = imin(edgecnt, BlockPerGrid);
		flip32_explore_new<<<blocks, (edgecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_edge, edgecnt, dev_halfedge, qualmeasure);
		cudaMemcpy(cedge, dev_edge, edgecnt*sizeof(cu_edge), cudaMemcpyDeviceToHost);

		// edge selecting
		selectedge = new int[edgecnt];
		edgeSelecting_new(cedge, edgecnt, tetcnt, selectedge, selectedgecnt);

		// flip32 (set -1 to the meshtets[] of invalid tets)
		if (selectedgecnt)
		{
			int *dev_selectedge;
			cudaMalloc((void**)&dev_selectedge, selectedgecnt*sizeof(int));
			cudaMemcpy(dev_selectedge, selectedge, selectedgecnt*sizeof(int), cudaMemcpyHostToDevice);

			// do flip32
			blocks = imin(selectedgecnt, BlockPerGrid);
			//flip32<<<blocks, (selectedgecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_edge, dev_selectedge, selectedgecnt, qualmeasure);

			flipsucc += selectedgecnt;
			cudaFree(dev_selectedge);
		}
		-- loop;
	}
}

/************************** end of flipping ******************************/


/************************ New Flip ****************************/
__device__ void add_vec3f(float *a, float *b, float *res)
{
	for (int i = 0; i < 3; i++)
		res[i] = a[i] + b[i];
}

__device__ void minus_vec3f(float *a, float *b, float *res)
{
	for (int i = 0; i < 3; i++)
		res[i] = a[i] - b[i];
}

__device__ float dot_vec3f(float *a, float *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

__device__ void cross_vec3f(float *a, float *b, float *res)
{
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ float norm_vec3f(float *a)
{
	return sqrtf(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

__device__ void normalize_vec3f(float *a)
{
	float norm = norm_vec3f(a);
	for (int i = 0; i < 3; i++)
		a[i] /= norm;
}

__device__ int checkPointInsideTriangle(float triangle[3][3], float p[3])
{
	float edge1[3];
	float edge2[3];
	float faceNormal[3];

	minus_vec3f(triangle[1], triangle[0], edge1);
	minus_vec3f(triangle[2], triangle[0], edge2);
	cross_vec3f(edge1, edge2, faceNormal);
	normalize_vec3f(faceNormal);

	float v[3][3];
	float vec1[3], vec2[3];

	minus_vec3f(triangle[1], triangle[0], vec1);
	minus_vec3f(p, triangle[0], vec2);
	cross_vec3f(vec1, vec2, v[0]);

	minus_vec3f(triangle[2], triangle[1], vec1);
	minus_vec3f(p, triangle[1], vec2);
	cross_vec3f(vec1, vec2, v[1]);

	minus_vec3f(triangle[0], triangle[2], vec1);
	minus_vec3f(p, triangle[2], vec2);
	cross_vec3f(vec1, vec2, v[2]);

	float cosAngle;

	for (int i = 0; i < 3; i ++)
	{
		if (norm_vec3f(v[i]) > 1e-10)
		{
			cosAngle = dot_vec3f(v[i], faceNormal) / norm_vec3f(v[i]);
		}
		else
		{
			cosAngle = 0;
		}
		if (cosAngle < - 1e-10)
		{
			// outside the triangle
			return 0;
		}
		else if (cosAngle == 0)
		{
			// on the triangle boundary
			return -1;
		}
	}
	// inside the triangle
	return 1;
}

__device__ float vertexTriangleDistance(float *points, int triangle[3], int p_, int & pStatus)
{
	float dis;
	dis = 0;

	float trianglev[3][3];
	float edge1[3], edge2[3];
	float faceNormal[3];
	float norm;
	float p[3];
	float pp[3];

	// get p
	for(int i = 0; i < 3; i ++)
		p[i] = points[p_*3+i];

	// get triangle points
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			trianglev[i][j] = points[triangle[i]*3+j];

	minus_vec3f(trianglev[1], trianglev[0], edge1);
	minus_vec3f(trianglev[2], trianglev[0], edge2);
	cross_vec3f(edge1, edge2, faceNormal);
	normalize_vec3f(faceNormal);
	norm = norm_vec3f(faceNormal);
	pStatus = 0;

	dis = (faceNormal[0] * p[0] + faceNormal[1] * p[1] +
		faceNormal[2] * p[2] - faceNormal[0] * trianglev[0][0] -
		faceNormal[1] * trianglev[0][1] - faceNormal[2] * trianglev[0][2]) / norm;

	for(int i = 0; i < 3; i ++)
		pp[i] = p[i] - faceNormal[i]*dis;

	pStatus = checkPointInsideTriangle(trianglev, pp);
	return dis;
}

__device__ bool trianglehasedge(int ph[3], int e[2])
{
	int c = 0;
	for (int i = 0; i < 3; i++)
	{
		if (ph[i] == e[0] || ph[i] == e[1])
			++ c;
	}
	if (c == 2)
		return true;
	return false;
}

__device__ bool trianglehaspoint(int ph[3], int p)
{
	return (ph[0] == p || ph[1] == p || ph[2] == p);
}

/* 确认风筝型四面体周围的四面体分布是否复合要求 */
__device__ bool kitesituationcheck(float *points, cu_tetra *tetra, int tidx, cu_flip23face *face, cu_halfface *halfface, 
								   int e1[2], int e2[2], int qualmeasure)
{
	int hf1[2];
	int hf2[2];
	int ophf1[2];
	int ophf2[2];
	cu_tetra tet1[2];
	cu_tetra tet2[2];
	int tetidx[4];
	int p1[2], p2[2];
	int fidx1, fidx2;
	cu_halfface tmphface;
	cu_flip23face tmpface;

	fidx1 = fidx2 = 0;
	for (int i = 0; i < 4; i ++)
	{
		tmphface = halfface[tidx*4+i];
		if (trianglehasedge(tmphface.pointhandle, e1))
			hf1[fidx1++] = tidx*4+i;
		else if(trianglehasedge(tmphface.pointhandle, e2))
			hf2[fidx2++] = tidx*4+i;
	}

	// find opposite halfface
	tmpface = face[halfface[hf1[0]].face];
	if (tmpface.hf[0] == hf1[0])
		ophf1[0] = tmpface.hf[1];
	else
		ophf1[0] = tmpface.hf[0];

	tmpface = face[halfface[hf1[1]].face];
	if (tmpface.hf[0] == hf1[1])
		ophf1[1] = tmpface.hf[1];
	else
		ophf1[1] = tmpface.hf[0];

	tmpface = face[halfface[hf2[0]].face];
	if (tmpface.hf[0] == hf2[0])
		ophf2[0] = tmpface.hf[1];
	else
		ophf2[0] = tmpface.hf[0];

	tmpface = face[halfface[hf2[1]].face];
	if (tmpface.hf[0] == hf2[1])
		ophf2[1] = tmpface.hf[1];
	else
		ophf2[1] = tmpface.hf[0];

	// get top points
	for (int i = 0; i < 2; i++)
	{
		tet1[i] = tetra[ophf1[i]>>2];
		tetidx[i] = ophf1[i]>>2;
		tmphface = halfface[ophf1[i]];
		for (int j = 0; j < 4; j++)
		{
			if (!trianglehaspoint(tmphface.pointhandle, tet1[i].v[j]))
			{
				p1[i] = tet1[i].v[j];
				break;
			}
		}

		tet2[i] = tetra[ophf2[i]>>2];
		tetidx[i+2] = ophf2[i]>>2;
		tmphface = halfface[ophf2[i]];
		for (int j = 0; j < 4; j++)
		{
			if (!trianglehaspoint(tmphface.pointhandle, tet2[i].v[j]))
			{
				p2[i] = tet2[i].v[j];
				break;
			}
		}
	}

	if (p1[0] == p1[1] && p2[0] == p2[1])
	{
		float p[4][3];
		cu_tetra currtet;
		int newtet[16];
		float qualtet1[2], qualtet2[2], tetqual;
		float qualnewtet[2];
		float minqualbefore, minqualafter1, minqualafter2;

		// calculate quality before newflip
		minqualbefore = 1.0;
		for (int k = 0; k < 2; k ++)
		{
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 3; j++)
					p[i][j] = points[tet1[k].v[i]*3+j];
			qualtet1[k] = tetquality(p, qualmeasure);
			if (minqualbefore > qualtet1[k])
				minqualbefore = qualtet1[k];

			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 3; j++)
					p[i][j] = points[tet2[k].v[i]*3+j];
			qualtet2[k] = tetquality(p, qualmeasure);
			if (minqualbefore > qualtet2[k])
				minqualbefore = qualtet2[k];
		}

		currtet = tetra[tidx];
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				p[i][j] = points[currtet.v[i]*3+j];

		tetqual = tetquality(p, qualmeasure);
		if (minqualbefore > tetqual)
			minqualbefore = tetqual;

		// update tetras' information
		tetra[tidx].quality = minqualbefore;
		tetra[tidx].fliptype = 1;
		tetra[tidx].newflipface[0] = hf1[0];
		tetra[tidx].newflipface[1] = hf1[1];
		tetra[tidx].newflipface[2] = hf2[0];
		tetra[tidx].newflipface[3] = hf2[1];
		tetra[tidx].tet[0] = tetidx[0];
		tetra[tidx].tet[1] = tetidx[1];
		tetra[tidx].tet[2] = tetidx[2];
		tetra[tidx].tet[3] = tetidx[3];
		tetra[tidx].flippoint[0] = p1[0];
		tetra[tidx].flippoint[1] = p2[0];

		// 确定flip策略
		// 组合策略一：e1周围两个四面体+新的两个四面体
		newtet[0] = p2[0];
		newtet[1] = halfface[hf1[0]].pointhandle[0];
		newtet[2] = halfface[hf1[0]].pointhandle[1];
		newtet[3] = halfface[hf1[0]].pointhandle[2];

		newtet[4] = p2[0];
		newtet[5] = halfface[hf1[1]].pointhandle[0];
		newtet[6] = halfface[hf1[1]].pointhandle[1];
		newtet[7] = halfface[hf1[1]].pointhandle[2];

		for (int k = 0; k < 2; k ++)
		{
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 3; j++)
					p[i][j] = points[newtet[k*4+i]*3+j];
			qualnewtet[k] = tetquality(p, qualmeasure);
		}

		minqualafter1 = 1.0;
		for (int i = 0; i < 2; i ++)
		{
			if (minqualafter1 > qualtet1[i])
				minqualafter1 = qualtet1[i];
			if (minqualafter1 > qualnewtet[i])
				minqualafter1 = qualnewtet[i];
		}

		// 组合策略二：e2周围两个四面体+新的两个四面体
		newtet[8]  = p1[0];
		newtet[9]  = halfface[hf2[0]].pointhandle[0];
		newtet[10] = halfface[hf2[0]].pointhandle[1];
		newtet[11] = halfface[hf2[0]].pointhandle[2];

		newtet[12] = p1[0];
		newtet[13] = halfface[hf2[1]].pointhandle[0];
		newtet[14] = halfface[hf2[1]].pointhandle[1];
		newtet[15] = halfface[hf2[1]].pointhandle[2];

		for (int k = 2; k < 4; k ++)
		{
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 3; j++)
					p[i][j] = points[newtet[k*4+i]*3+j];
			qualnewtet[k] = tetquality(p, qualmeasure);
		}

		minqualafter2 = 1.0;
		for (int i = 0; i < 2; i ++)
		{
			if (minqualafter2 > qualtet2[i])
				minqualafter2 = qualtet2[i];
			if (minqualafter2 > qualnewtet[i])
				minqualafter2 = qualnewtet[i];
		}

		// update information
		if (minqualafter2 > minqualafter1 && minqualafter2 > minqualbefore)
		{
			tetra[tidx].strategy = 2;
			tetra[tidx].val = minqualafter2;
		}
		else if (minqualafter1 > minqualafter2 && minqualafter1 > minqualbefore)
		{
			tetra[tidx].strategy = 1;
			tetra[tidx].val = minqualafter1;
		}
		return true;
	}
	return false;
}

/* 确认三角型四面体周围的四面体分布是否复合要求 */
__device__ bool trianglesituationcheck(float *points, cu_tetra *tetra, int tidx, cu_flip23face *face, cu_halfface *halfface, int vidx, int qualmeasure)
{
	// 与vh_周围三个面相邻的三个四面体相交于一点
	int hf[3];
	hf[0] = (tidx<<2)+(vidx+1)%4;
	hf[1] = (tidx<<2)+(vidx+2)%4;
	hf[2] = (tidx<<2)+(vidx+3)%4;

	int ophf[3];
	int ph[3];
	int tetidx[3];
	cu_flip23face tmpface;
	cu_tetra tet[3];
	cu_tetra currtet;
	for (int i = 0; i < 3; i ++)
	{
		// get opposite halfface
		tmpface = face[halfface[hf[i]].face];
		if (tmpface.hf[0] == hf[i])
			ophf[i] = tmpface.hf[1];
		else
			ophf[i] = tmpface.hf[0];

		// get the tetra contains the opposite halfface
		tet[i] = tetra[ophf[i]>>2];
		tetidx[i] = ophf[i]>>2;
		// get the top point
		ph[i] = tet[i].v[ophf[i]%4]; 
	}

	if (ph[0] == ph[1] && ph[0] == ph[2])
	{
		float p[4][3];
		float minqualbefore, minqualafter, tetqual;
		int newtet[4];
		cu_halfface tmphf;
		int hface;

		// calculate quality before newflip
		minqualbefore = 1.0;
		for (int k = 0; k < 3; k ++)
		{
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 3; j++)
					p[i][j] = points[tet[k].v[i]*3+j];
			tetqual = tetquality(p, qualmeasure);
			if (minqualbefore > tetqual)
				minqualbefore = tetqual;
		}

		currtet = tetra[tidx];
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				p[i][j] = points[currtet.v[i]*3+j];
		tetqual = tetquality(p, qualmeasure);
		if (minqualbefore > tetqual)
			minqualbefore = tetqual;

		tetra[tidx].quality = minqualbefore;
		tetra[tidx].fliptype = 2;
		tetra[tidx].newflipface[0] = (tidx<<2)+vidx;
		tetra[tidx].newflipface[1] = hf[0];
		tetra[tidx].newflipface[2] = hf[1];
		tetra[tidx].newflipface[3] = hf[2];
		tetra[tidx].tet[0] = tetidx[0];
		tetra[tidx].tet[1] = tetidx[1];
		tetra[tidx].tet[2] = tetidx[2];
		tetra[tidx].flippoint[0] = ph[0];

		// newtet
		hface = (tidx<<2)+vidx;
		tmphf = halfface[hface];
		newtet[0] = ph[0];
		newtet[1] = tmphf.pointhandle[0];
		newtet[2] = tmphf.pointhandle[1];
		newtet[3] = tmphf.pointhandle[2];

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				p[i][j] = points[newtet[i]*3+j];
		minqualafter = tetquality(p, qualmeasure);

		if (minqualafter > minqualbefore)
			tetra[tidx].val = minqualafter;
		return true;
	}
	return false;
}

__global__ void newflip_explore(float *points, int pointcnt, cu_tetra *tetra, int tetcnt, cu_flip23face *face, int facecnt, 
					            cu_halfface *halfface, int halffacecnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;

	int p[4];
	int trianglepoint[3];
	int pp;
	int flag;
	int state[3];         // 顶点位置统计 0: outside   1: inside   2: on the boundary
	int vhandle;             // 三角形情况的顶点
	int situation;        // 0：风筝型 1：三角形  -1: 其他
	int e1[2], e2[2];
	//bool issuitable;
	cu_tetra currtet;
	//int count = 0;

	while(tid < tetcnt)
	{
		currtet = tetra[tid];
		if (currtet.isboundary)
		{
			tid += offset;
			continue;
		}

		// get local data
		for (int i = 0; i < 4; i++)
			p[i] = currtet.v[i];

		situation = 0;
		state[0] = state[1] = state[2] = 0;

		/* 确定四面体形状 */
		for (int i = 0; i < 4; i++)
		{
			// 将一点投影到另外三个点所在的平面
			pp = p[i];

			for (int j = 0; j < 3; j++)
				trianglepoint[j] = p[(i+j+1)%4];

			vertexTriangleDistance(points, trianglepoint, pp, flag);

			if (flag == 0)
				++state[0];
			else if (flag == 1)
			{
				++state[1];
				vhandle = i;    // 顶点编号记录
			}
			else
				++state[2];
		}

		if (state[1] == 0)
			situation = 0;
		else if (state[1] == 1)
			situation = 1;
		else
			situation = -1;

		// 确定四面体周围四面体分布情况
		if(situation == 0)//风筝型
		{
// 			cu_tetra tet;
// 			tet = tetra[tid];

			// get tet points
			int v[4];
			//float p[4][3];
			//for (int i = 0; i < 4; i++)
			//{
			//	v[i] = tetra[tid].v[i];
			//	for (int j = 0; j < 3; j++)
			//		p[i][j] = points[v[i]*3+j];
			//}

			// get tet edges
			int edge[6][2];
			edge[0][0] = 0; edge[0][1] = 1;
			edge[1][0] = 0; edge[1][1] = 2;
			edge[2][0] = 0; edge[2][1] = 3;
			edge[3][0] = 2; edge[3][1] = 3;
			edge[4][0] = 1; edge[4][1] = 3;
			edge[5][0] = 1; edge[5][1] = 2;

			for (int i = 0; i < 3; i ++)
			{
				e1[0] = v[edge[i][0]];
				e1[1] = v[edge[i][1]];

				e2[0] = v[edge[(i+3)%6][0]];
				e2[1] = v[edge[(i+3)%6][1]];
				kitesituationcheck(points,tetra, tid, face, halfface, e1, e2, qualmeasure);
			}
		}
		else if(situation == 1)//三角型
		{
			trianglesituationcheck(points, tetra, tid, face, halfface, vhandle, qualmeasure);
		}
		tid += offset;
	}
}

void newfliptetraSelecting(struct cu_tetra *tetra, int tetcnt, int *selecttet, int &selecttetcnt)
{
	int stetcnt = 0;
	int* stet;
	int i, j;
	//float threashold = 1.0e-5;

	stet = new int[tetcnt];
	/* pick out the tets by flipping succeed val and sort it*/
	for (i = 0; i < tetcnt; i++)
	{
		if (tetra[i].val < 0 || !(tetra[i].val > tetra[i].quality))
			continue;

		// insert the tet
		for (j = stetcnt-1; j > -1 ; j--)
		{
			if (tetra[stet[j]].quality > tetra[i].quality)
				stet[j+1] = stet[j];
			else 
				break;
		}
		stet[j+1] = i;
		++ stetcnt;
	}

	/* if there are no any face meeting the requirement */
	if (!stetcnt)
	{
		selecttetcnt = 0;
		selecttet = NULL;
		return;
	}

	/* select a new tet set */
	int tet[4];
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, 1, tetcnt*sizeof(bool));

	// push the first tet
	selecttet[0] = stet[0];
	selecttetcnt = 1;

	// set flags of tets incident to the tet
	tetflag[selecttet[0]] = 0;
	if (tetra[selecttet[0]].fliptype == 1)
	{
		for (i = 0; i < 4; i ++)
			tetflag[tetra[selecttet[0]].tet[i]] = 0;
	}
	else
	{
		for (i = 0; i < 3; i ++)
			tetflag[tetra[selecttet[0]].tet[i]] = 0;
	}

	/* if the tets incident to a face are available, 
	   then add the face into the array*/
	for (i = 1; i < stetcnt; i++)
	{
		if (tetra[stet[i]].fliptype == 1)
		{
			for (j = 0; j < 4; j ++)
				tet[j] = tetra[stet[i]].tet[j];
			if (tetflag[tet[0]] && tetflag[tet[1]] && tetflag[tet[2]] && tetflag[tet[3]] && tetflag[stet[i]])
			{
				selecttet[selecttetcnt++] = stet[i];
				tetflag[stet[i]] = 0;
				tetflag[tet[0]] = 0;
				tetflag[tet[1]] = 0;
				tetflag[tet[2]] = 0;
				tetflag[tet[3]] = 0;
			}
		}
		else
		{
			for (j = 0; j < 3; j ++)
				tet[j] = tetra[stet[i]].tet[j];
			if (tetflag[tet[0]] && tetflag[tet[1]] && tetflag[tet[2]] && tetflag[stet[i]])
			{
				selecttet[selecttetcnt++] = stet[i];
				tetflag[stet[i]] = 0;
				tetflag[tet[0]] = 0;
				tetflag[tet[1]] = 0;
				tetflag[tet[2]] = 0;
			}
		}
	}
}

__global__ void newflip(cu_tetra *tetra, int *selecttet, int selecttetcnt, cu_halfface *halfface)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;

	int fliptype;
	cu_tetra currtet;
	//float p[4][3];
	//float qual;
	//float newqual[2];
	while(tid < selecttetcnt)
	{
		currtet = tetra[selecttet[tid]];
		fliptype = currtet.fliptype;

		// 风筝型
		if (fliptype == 1)
		{
			int hf1[2], hf2[2];
			int tet1[2], tet2[2];
			int p1, p2;
			int newtet[8];

			hf1[0] = currtet.newflipface[0];
			hf1[1] = currtet.newflipface[1];
			hf2[0] = currtet.newflipface[2];
			hf2[1] = currtet.newflipface[3];
			tet1[0] = currtet.tet[0];
			tet1[1] = currtet.tet[1];
			tet2[0] = currtet.tet[2];
			tet2[1] = currtet.tet[3];
			p1 = currtet.flippoint[0];
			p2 = currtet.flippoint[1];

			if (currtet.strategy == 1)
			{
				newtet[0] = p2;
				newtet[1] = halfface[hf1[0]].pointhandle[0];
				newtet[2] = halfface[hf1[0]].pointhandle[1];
				newtet[3] = halfface[hf1[0]].pointhandle[2];

				newtet[4] = p2;
				newtet[5] = halfface[hf1[1]].pointhandle[0];
				newtet[6] = halfface[hf1[1]].pointhandle[1];
				newtet[7] = halfface[hf1[1]].pointhandle[2];

				for (int i = 0; i < 2; i ++)
				{
					for (int j = 0; j < 4; j++)
						tetra[tet2[i]].v[j] = newtet[i*4+j];
				}
			}
			else
			{
				newtet[0] = p1;
				newtet[1] = halfface[hf2[0]].pointhandle[0];
				newtet[2] = halfface[hf2[0]].pointhandle[1];
				newtet[3] = halfface[hf2[0]].pointhandle[2];

				newtet[4] = p1;
				newtet[5] = halfface[hf2[1]].pointhandle[0];
				newtet[6] = halfface[hf2[1]].pointhandle[1];
				newtet[7] = halfface[hf2[1]].pointhandle[2];

				for (int i = 0; i < 2; i ++)
				{
					for (int j = 0; j < 4; j++)
						tetra[tet1[i]].v[j] = newtet[i*4+j];
				}
			}
		}
		// 三角型
		else
		{
			int ph;
			int hf;
			int tetidx[3];
			int newtet[4];

			ph = currtet.flippoint[0];
			hf = currtet.newflipface[0];
			for (int i = 0; i < 3; i++)
				tetidx[i] = currtet.tet[i];

			// newtet
			newtet[0] = ph;
			newtet[1] = halfface[hf].pointhandle[0];
			newtet[2] = halfface[hf].pointhandle[1];
			newtet[3] = halfface[hf].pointhandle[2];

			for (int i = 0; i < 4; i++)
			{
				tetra[tetidx[0]].v[i] = newtet[i];
				tetra[tetidx[1]].v[i] = -1;
				tetra[tetidx[2]].v[i] = -1;
			}
		}
		for (int i = 0; i < 4; i++)
			tetra[tid].v[i] = -1;

		tid += offset;
	}
}

extern "C" void cuda_newflip(float *points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
							 int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
							 int &flipsucc, float &time)
{
	// create face and halfface structure
	struct cu_flip23face *cface;
	struct cu_halfface *chalfface;
	struct cu_tetra *ctetra;
	float qualbefore;
	float qualafter;
	//float qualtmp;
	int tetcapacity;
	int facecapacity;
	int halffacecapacity;

	tetcapacity = int(tetcnt);
	facecapacity = int(facecnt);
	halffacecapacity = int(halffacecnt);
	ctetra = new struct cu_tetra[tetcapacity];
	cface = new struct cu_flip23face[facecapacity];
	chalfface = new struct cu_halfface[halffacecapacity];

	// get face
	for (int i = 0; i < facecnt; i++)
	{
		cface[i].hf[0] = face[2*i];
		cface[i].hf[1] = face[2*i+1];
		cface[i].quality = -1;
		cface[i].val = -1;
	}

	// get halfface
	for (int i = 0; i < halffacecnt; i++)
	{
		chalfface[i].pointhandle[0] = halfface[4*i];
		chalfface[i].pointhandle[1] = halfface[4*i+1];
		chalfface[i].pointhandle[2] = halfface[4*i+2];
		chalfface[i].face = halfface[4*i+3];
	}

	// get tetra
	for (int i = 0; i < tetcnt; i ++)
	{
		ctetra[i].v[0] = meshtets[5*i];
		ctetra[i].v[1] = meshtets[5*i+1];
		ctetra[i].v[2] = meshtets[5*i+2];
		ctetra[i].v[3] = meshtets[5*i+3];
		ctetra[i].isboundary = meshtets[5*i+4];
		ctetra[i].fliptype = -1;
		ctetra[i].newflipvertex = -1;
		ctetra[i].quality = -1;
		ctetra[i].val = -1;
		for (int j = 0; j < 4; j ++)
			ctetra[i].newflipface[j] = 0;
	}

	// calculate quality before flipping
	cutetramesh_quality(points, pointcnt, ctetra, tetcnt, qualmeasure, qualbefore);

	// CUDA Parallel
	// 分配设备存储空间
	//int loop = 1;
	float *dev_points;
	//int *dev_meshtets, *dev_tempmeshtets;
	struct cu_flip23face *dev_face;
	struct cu_halfface *dev_halfface;
	struct cu_tetra *dev_tetra;
	int *selecttetra;
	int selecttetracnt;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// flip23 会生成新的四面体、半面和面，因此先预分配一部分空间
	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
	//cudaMalloc((void**)&dev_meshtets, 4*tetcapacity*sizeof(int));
	cudaMalloc((void**)&dev_face, facecapacity*sizeof(cu_flip23face));
	cudaMalloc((void**)&dev_halfface, halffacecapacity*sizeof(cu_halfface));
	cudaMalloc((void**)&dev_tetra, tetcapacity*sizeof(cu_tetra));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_face, cface, facecnt*sizeof(cu_flip23face), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_halfface, chalfface, halffacecnt*sizeof(cu_halfface), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_tetra, ctetra, tetcnt*sizeof(cu_tetra), cudaMemcpyHostToDevice);

	// do new flip
	// newflip explore
	int blocks = imin(tetcnt, BlockPerGrid);
	newflip_explore<<<blocks, (tetcnt+blocks-1)/blocks>>>(dev_points, pointcnt, dev_tetra, tetcnt, dev_face, facecnt, 
		dev_halfface, halffacecnt, qualmeasure); 
	cudaMemcpy(ctetra, dev_tetra, tetcnt*sizeof(cu_tetra), cudaMemcpyDeviceToHost);

	// tetra selecting
	selecttetra = new int[tetcnt];
	selecttetracnt = 0;
	newfliptetraSelecting(ctetra, tetcnt, selecttetra, selecttetracnt);

	// newflip (set -1 to the meshtets[] of invalid tets)
	if (selecttetracnt)
	{
		int *dev_selecttetra;
		cudaMalloc((void**)&dev_selecttetra, selecttetracnt*sizeof(int));
		cudaMemcpy(dev_selecttetra, selecttetra, selecttetracnt*sizeof(int), cudaMemcpyHostToDevice);

		// do newflip
		blocks = imin(selecttetracnt, BlockPerGrid);
		newflip<<<blocks, (selecttetracnt+blocks-1)/blocks>>>(dev_tetra, dev_selecttetra, selecttetracnt, dev_halfface);

		flipsucc += selecttetracnt;
		cudaFree(dev_selecttetra);
	}		

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);

	time = elaspsedTime;

	// calculate quality after flipping
	qualafter = 1.0;
	cudaMemcpy(ctetra, dev_tetra, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToHost);
	// calculate quality before flipping
	cutetramesh_quality(points, pointcnt, ctetra, tetcnt, qualmeasure, qualafter);

	qualbefore_ = qualbefore;
	qualafter_ = qualafter;
	flipsucc = selecttetracnt;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_points);
	//cudaFree(dev_meshtets);
	cudaFree(dev_halfface);
	cudaFree(dev_face);

	delete [] cface;
	delete [] chalfface;
	delete [] ctetra;
	delete [] selecttetra;
}
/************************ End of New Flip *****************************/

/************************ Edge Contraction ****************************/

__device__ void getECNewTetras(int *newtetra, int *newtetracnt, int *fromtetra, int fromtetracnt, int *totetra, int totetracnt,
							   int *edgestar, int edgestarcnt)
{
	*newtetracnt = 0;
	int i, j, k;
	int currft, currtt;
	for (i = 0, j = 0; i < fromtetracnt || j < totetracnt; i++, j++)
	{
		if (i < fromtetracnt)
			currft = fromtetra[i];
		else
			currft = -1;

		if (j < totetracnt)
			currtt = totetra[j];
		else
			currtt = -1;

		for (k = 0; k < edgestarcnt; k++)
		{
			if (edgestar[k] == fromtetra[i])
				currft = -1;
			if(edgestar[k] == totetra[j])
				currtt = -1;
		}

		if (currft != -1)
			newtetra[(*newtetracnt)++] = currft;
		if (currtt != -1)
			newtetra[(*newtetracnt)++] = currtt;
	}
}

// for model: out, rand1, rand2
//#define MAXPOINTSTAR 300
//#define MAXPOINTSTARSUM 550

// for other models
#define MAXPOINTSTAR 70
#define MAXPOINTSTARSUM 120

#define MASEDGESTAR 25

// not tested yet
__global__ void edge_contraction_explore(float *points, int *meshtets, cu_edge *edge, int edgecnt, cu_halfedge *halfedge, 
										 int *incidenttet, int *incidenttetcnt, int largesttetcnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	int fromv, tov, tempv;
	float fp[3], tp[3], mp[3];
	int edgestar[MASEDGESTAR];
	int newtetra[MAXPOINTSTARSUM];
	int fromtetra[MAXPOINTSTAR];
	int totetra[MAXPOINTSTAR];
	int edgestarcnt;
	int newtetracnt;
	int fromtetracnt;
	int totetracnt;
	float newtetrapoint[MAXPOINTSTARSUM*12];
	int *phe;
	float tempqual;
	float qualbefore, qualafter;
	cu_edge curredge;
	int i, j, k;

	while(tid < edgecnt)
	{
		curredge = edge[tid];

		// if the edge is boundary go next
		if (curredge.is_boundary) 
		{
			tid += offset;
			continue;
		}

		/* fetch local data */
		// get edge star
		phe = curredge.halfedge;
		edgestarcnt = curredge.halfedgecnt / 2;
		for (i = 0; i < edgestarcnt; i ++)
			edgestar[i] = (phe[i*2])/12;

		// get new tetras
		fromv = halfedge[phe[0]].fromv;
		tov = halfedge[phe[0]].tov;
		fromtetracnt = incidenttetcnt[fromv];
		totetracnt = incidenttetcnt[tov];

		for (i = 0; i < fromtetracnt; i++)
			fromtetra[i] = incidenttet[fromv*largesttetcnt + i];
		for (i = 0; i < totetracnt; i++)
			totetra[i] = incidenttet[tov*largesttetcnt + i];

		getECNewTetras(newtetra, &newtetracnt, fromtetra, fromtetracnt, totetra, totetracnt, edgestar, edgestarcnt);

		// calculate quality
		qualbefore = minstackquality(points, meshtets, newtetra, newtetracnt, qualmeasure);
		tempqual = minstackquality(points, meshtets, edgestar, edgestarcnt, qualmeasure);
		qualbefore = qualbefore < tempqual ? qualbefore : tempqual;

		// get endpoints
		for (i = 0; i < 3; i ++)
		{
			fp[i] = points[3*fromv+i];
			tp[i] = points[3*tov + i];
			mp[i] = (fp[i] + tp[i])/2.0;
		}

		// get newtetras' points
		for (i = 0; i < newtetracnt; i++)
		{
			for (j = 0; j < 4; j++)
			{
				tempv = meshtets[4*newtetra[i]+j];
				if (tempv == fromv || tempv == tov)
				{
					for (k = 0; k < 3; k++)
						newtetrapoint[3*(i*4+j)+k] = mp[k];
				}
				else
				{
					for (k = 0; k < 3; k++)
						newtetrapoint[3*(i*4+j)+k] = points[3*tempv+k];
				}
			}
		}

		// calculate new quality
		qualafter = minstackquality(newtetrapoint, newtetracnt, qualmeasure);

		// record data
		edge[tid].quality = qualbefore;
		edge[tid].val = qualafter - qualbefore;

		tid += offset;
	}
}


__global__ void edge_contraction(float *points, int *meshtets, cu_edge *edge, cu_halfedge *halfedge,
								 int *selectedge, int selectedgecnt)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;

	int fromv, tov;
	float fp[3], tp[3], mp[3];
	cu_edge curredge;
	int tetidx;
	int i, j;

	while(tid < selectedgecnt)
	{
		curredge = edge[selectedge[tid]];

		// set two endpoints with new coordinates
		fromv = halfedge[curredge.halfedge[0]].fromv;
		tov = halfedge[curredge.halfedge[0]].tov;
		for (i = 0; i < 3; i++)
		{
			fp[i] = points[3*fromv+i];
			tp[i] = points[3*tov+i];
			mp[i] = (fp[i] + tp[i])/2.0;
			points[3*fromv+i] = mp[i];
			points[3*tov+i] = mp[i];
		}

		// set invalid data to tetras incident to current edge
		for (i = 0; i < curredge.halfedgecnt; i++)
		{
			tetidx = (curredge.halfedge[i])/4;
			for (j = 0; j < 4; j++)
			{
				meshtets[4*tetidx+j] = -1;
			}
		}
		tid += offset;
	}
}

void edgeSelecting_EC(cu_edge *cuedge, int edgecnt, int tetcnt, int *selectedge, int &selectedgecnt)
{
	int sedgecnt = 0;
	int* sedge = new int[edgecnt];
	int i, j;
	float threashold = 1.0e-5;
	bool isavailable;
	cu_edge curredge;

	/* pick out the edges by flipping succeed val
	and sort it*/
	for (i = 0; i < edgecnt; i++)
	{
		if (!(cuedge[i].val > threashold))
			continue;

		// insert the edge
		for (j = sedgecnt-1; j > -1 ; j--)
		{
			if (cuedge[sedge[j]].quality > cuedge[i].quality)
				sedge[j+1] = sedge[j];
			else 
				break;
		}
		sedge[j+1] = i;
		++ sedgecnt;
	}

	/* if there are no any edge meeting the requirement */
	if (!sedgecnt)
	{
		selectedgecnt = 0;
		selectedge = NULL;
		return;
	}

	/* select a new edge set in which any two of them not in a same tet */
	//int tet;
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, true, tetcnt*sizeof(bool));

	// push the first edge
	selectedge[0] = sedge[0];
	selectedgecnt = 1;

	// set flags of tets incident to the first edge
	curredge = cuedge[sedge[0]];
	for (i = 0; i < curredge.halfedgecnt; i++)
	{
		tetflag[(curredge.halfedge[i])/12] = false;
	}

	/* if the tets incident to a edge are available, 
	then add the edge into the array*/
	for (i = 1; i < sedgecnt; i++)
	{
		isavailable = true;
		curredge = cuedge[sedge[i]];
		for (j = 0; j < curredge.halfedgecnt; j++)
		{
			if (!tetflag[(curredge.halfedge[j])/12])
			{
				isavailable = false;
				break;
			}
		}
		if (isavailable)
		{
			selectedge[selectedgecnt++] = sedge[i];
			for (j = 0; j < curredge.halfedgecnt; j++)
			{
				tetflag[(curredge.halfedge[j])/12] = false;
			}
		}
	}
	delete [] tetflag;
}

void edgeSelecting_EC_1(cu_edge *cuedge, int edgecnt, cu_halfedge *cuhalfedge, int halfedgecnt, int tetcnt, 
					  int *incidenttet, int *incidenttetcnt, int largesttetcnt, int *selectedge, int &selectedgecnt)
{
	int sedgecnt = 0;
	int* sedge = new int[edgecnt];
	int i, j;
	int fromv, tov;
	float threashold = 1.0e-5;
	bool isavailable;
	cu_edge curredge;

	/* pick out the edges by flipping succeed val
	and sort it*/
	for (i = 0; i < edgecnt; i++)
	{
		if (!(cuedge[i].val > threashold))
			continue;

		// insert the edge
		for (j = sedgecnt-1; j > -1 ; j--)
		{
			if (cuedge[sedge[j]].quality > cuedge[i].quality)
				sedge[j+1] = sedge[j];
			else 
				break;
		}
		sedge[j+1] = i;
		++ sedgecnt;
	}

	/* if there are no any edge meeting the requirement */
	if (!sedgecnt)
	{
		selectedgecnt = 0;
		selectedge = NULL;
		return;
	}

	/* select a new edge set in which any two of them not in a same tet */
	//int tet;
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, true, tetcnt*sizeof(bool));

	// push the first edge
	selectedge[0] = sedge[0];
	selectedgecnt = 1;

	// set flags of tets incident to the first edge
	curredge = cuedge[sedge[0]];
	fromv = cuhalfedge[curredge.halfedge[0]].fromv;
	tov = cuhalfedge[curredge.halfedge[0]].tov;
	for (i = 0; i < incidenttetcnt[fromv]; i++)
	{
		tetflag[incidenttet[fromv*largesttetcnt+i]] = false;
	}
	for (i = 0; i < incidenttetcnt[tov]; i++)
	{
		tetflag[incidenttet[tov*largesttetcnt+i]] = false;
	}

	/* if the tets incident to a edge are available, 
	then add the edge into the array*/
	for (i = 1; i < sedgecnt; i++)
	{
		isavailable = true;
		curredge = cuedge[sedge[i]];
		fromv = cuhalfedge[curredge.halfedge[0]].fromv;
		tov = cuhalfedge[curredge.halfedge[0]].tov;
		for (j = 0; j < incidenttetcnt[fromv] && isavailable; j++)
		{
			if (!tetflag[incidenttet[fromv*largesttetcnt+j]])
				isavailable = false;
		}
		for (j = 0; j < incidenttetcnt[tov] && isavailable; j++)
		{
			if (!tetflag[incidenttet[tov*largesttetcnt+j]])
				isavailable = false;
		}

		if (isavailable)
		{
			selectedge[selectedgecnt++] = sedge[i];
			for (j = 0; j < incidenttetcnt[fromv]; j++)
			{
				tetflag[incidenttet[fromv*largesttetcnt+j]] = false;
			}
			for (j = 0; j < incidenttetcnt[tov]; j++)
			{
				tetflag[incidenttet[tov*largesttetcnt+j]] = false;
			}
		}
	}
	delete [] tetflag;
}


extern "C" void cuda_edgeContraction(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
									 int *halfedge, int halfedgecnt, int *incidenttet, int* incidenttetcnt, int &largesttetcnt, 
									 int qualmeasure, float &qualbefore_, float &qualafter_, int &succ_, float &time)
{
	//FILE *outfile = NULL;
	//outfile = fopen("F:\\kernel_data_output.txt", "w");

	cu_edge *cedge;
	cu_halfedge *chalfedge;
	cedge = new cu_edge[edgecnt];
	chalfedge = new cu_halfedge[halfedgecnt];


	float qualbefore, qualafter;
	float *dev_points;
	int *dev_meshtets;
	int *dev_incidenttet;
	int *dev_incidenttetcnt;
	struct cu_edge *dev_edge;
	struct cu_halfedge *dev_halfedge;
	int *selectedge;
	int selectedgecnt;
	int loop = 1;

	// calculate quality before edge contraction
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// build edge
	int idx = 0;
	for (int i = 0; i < edgecnt; i ++)
	{
		cedge[i].is_boundary = edge[idx++];
		cedge[i].halfedgecnt = edge[idx++];
		for (int j = 0; j < cedge[i].halfedgecnt; j++)
		{
			cedge[i].halfedge[j] = edge[idx++];
		}
		cedge[i].quality = cedge[i].val = -1;
	}

	// build halfedge
	for (int i = 0; i < halfedgecnt; i++)
	{
		chalfedge[i].edgehandle = halfedge[i*3];
		chalfedge[i].fromv = halfedge[i*3+1];
		chalfedge[i].tov = halfedge[i*3+2];
	}

 	cudaMalloc((void**)&dev_points, 3*pointcnt*sizeof(float));
 	cudaMalloc((void**)&dev_meshtets, 4*tetcnt*sizeof(int));
 	cudaMalloc((void**)&dev_edge, edgecnt*sizeof(cu_edge));
 	cudaMalloc((void**)&dev_halfedge, halfedgecnt*sizeof(cu_halfedge));
 	cudaMalloc((void**)&dev_incidenttet, pointcnt * largesttetcnt * sizeof(int));
 	cudaMalloc((void**)&dev_incidenttetcnt, pointcnt * sizeof(int));
 
 	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
 	cudaMemcpy(dev_meshtets, meshtets, 4*tetcnt*sizeof(int), cudaMemcpyHostToDevice);
 	cudaMemcpy(dev_edge, cedge, edgecnt*sizeof(cu_edge), cudaMemcpyHostToDevice);
 	cudaMemcpy(dev_halfedge, chalfedge, halfedgecnt*sizeof(cu_halfedge), cudaMemcpyHostToDevice);
 	cudaMemcpy(dev_incidenttet, incidenttet, pointcnt * largesttetcnt * sizeof(int), cudaMemcpyHostToDevice);
 	cudaMemcpy(dev_incidenttetcnt, incidenttetcnt, pointcnt * sizeof(int), cudaMemcpyHostToDevice);

	int succ = 0;

	//InitGPUSet();
	//cuPrintInit();
	while(loop)
	{
		/***********do parallel edge contraction**********/
		// edge contraction explore
		int blocks = imin(edgecnt, BlockPerGrid);
		edge_contraction_explore<<<blocks, (edgecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_edge, edgecnt, dev_halfedge, dev_incidenttet, dev_incidenttetcnt, largesttetcnt, qualmeasure);
		//cudaPrintfDisplay(outfile, true);
		//cudaPrintfEnd(); 
		//fclose(outfile);

		cudaMemcpy(cedge, dev_edge, edgecnt*sizeof(cu_edge), cudaMemcpyDeviceToHost);

		// edge selecting
		selectedge = new int[edgecnt];
		edgeSelecting_EC(cedge, edgecnt, tetcnt, selectedge, selectedgecnt);
		//edgeSelecting_EC_1(cedge, edgecnt, chalfedge, halfedgecnt, tetcnt, incidenttet, incidenttetcnt, largesttetcnt, selectedge, selectedgecnt);

		// edge contraction (set -1 to the meshtets[] of invalid tets)
		if (selectedgecnt)
		{
			int *dev_selectedge;
			cudaMalloc((void**)&dev_selectedge, selectedgecnt*sizeof(int));
			cudaMemcpy(dev_selectedge, selectedge, selectedgecnt*sizeof(int), cudaMemcpyHostToDevice);

			// do edge contraction
			blocks = imin(selectedgecnt, BlockPerGrid);
			edge_contraction<<<blocks, (selectedgecnt+blocks-1)/blocks>>>(dev_points, dev_meshtets, dev_edge, dev_halfedge, dev_selectedge, selectedgecnt);

			succ += selectedgecnt;
			cudaFree(dev_selectedge);
		}
		-- loop;
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);

	time = elaspsedTime;

	// calculate quality after flipping
	qualafter = 1.0;
	cudaMemcpy(meshtets, dev_meshtets, 4*tetcnt*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(points, dev_points, 3*pointcnt*sizeof(float), cudaMemcpyDeviceToHost);
	cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualafter);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_points);
	cudaFree(dev_meshtets);
	cudaFree(dev_halfedge);
	cudaFree(dev_edge);
	cudaFree(dev_incidenttet);
	cudaFree(dev_incidenttetcnt);

	delete [] cedge;
	delete [] chalfedge;

	qualbefore_ = qualbefore;
	qualafter_ = qualafter;
	succ_ = succ;
}
/************************ End of Edge Contraction ************************/

/************************ Parallel Vertex Insertion ***************************/

/* find tet adjacencies to a specified face 
* param[hf]: handle of the specified half face
* param[outin]: return the two adjacencies tets 
outin[0]: tet contains the hf's opposite face
outin[1]: tet contains the face hf
*/
__device__ int opposite_halfface(cu_face *face, cu_halfface *halfface, int hf)
{
	int fidx = halfface[hf].face;
	int oppf;
	if (face[fidx].hf[0] == hf)
		oppf = face[fidx].hf[1];
	else
		oppf = face[fidx].hf[0];
	return oppf;
}


/* find cavitytet by a tetrahedron handle*/
__device__ int findCavityTet(CavityTet *outcavity, int outcavitysize, int tethandle)
{
	int idx;
	for (idx = 0; idx < outcavitysize; idx++)
	{
		if (outcavity[idx].handle == tethandle) break;
	}
	return idx;
}

__device__ int find_array_int(int *array_, int arrsize, int data_)
{
	int i = 0;
	for (; i < arrsize; i++)
	{
		if (array_[i] == data_)
			break;
	}
	return i;
}

__device__ void halfFacePoints(float *point, cu_halfface *halfface, int hf, float fpoint[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		fpoint[i][0] = point[(halfface[hf].pointhandle[i])*3];
		fpoint[i][1] = point[(halfface[hf].pointhandle[i])*3+1];
		fpoint[i][2] = point[(halfface[hf].pointhandle[i])*3+2];
	}
}

/* check if a facet (indicated by fp1, fp2, fp3) is oriented to a vertex v */
__device__ int orient(float v[3], float fp1[3], float fp2[3], float fp3[3])
{
	float center[3];
	float facenormal[3];
	float result;
	float vec1[3], vec2[3];
	
	// get face center
	vector_add(fp1, fp2, center);
	vector_add(center, fp3, center);
	center[0] /= 3.0;
	center[1] /= 3.0;
	center[2] /= 3.0;

	vector_minus(fp2, fp1, vec1);
	vector_minus(fp3, fp1, vec2);
	vector_cross(vec1, vec2, facenormal);

	vector_minus(v, center, vec1);
	result = vector_dot(vec1, facenormal);
	if (result > 0)
		return 1;
	if (result == 0)
		return 0;
	return -1;
}

/* given a vertex, it's position, an initial cavitiy of tets and 
* a set of outward-oriented faces of that cavity, build a DAG
* representing the largest star-shaped cavity from the point of
* view of the inserted vertex */
__device__ void buildcavitydag(float *point, cu_InsertTet *tets, cu_face *face, cu_halfface *halfface, 
					           float pnew[3], int inserttetidx, cu_InsertTet &inserttet, float tetpoint[4][3], int qualmeasure, 
				               CavityTet *outcavity, int &outcavitysize)
{
	int F[MAXSTACKFACE];    /* candidate face list */
	int W[MAXSTACKFACE];    /* wall face list */
	int B[MAXSTACKTET];     /* blocking tet list */
	int Fsize=0, Wsize=0, Bsize=0;
	int currtet;                   /* current tet */
	int currf;                     /* current face */
	int oppf;
	int outin[2];
	int Wcount = 0;
	bool facing = true;
	int otherfaces[3];
	int deepest = 0;
	int nonwall[3];
	int numnonwall = 0;
	float fpoint[3][3];

	/* output cavity stuff */
	CavityTet cavtet;
	CavityFace cavface;
	CavityFace cavface2;
	int tetidx;
	int i,j;

	/* initialize cavity tet*/
	cavtet.handle = inserttetidx;
	cavtet.quality = tets[inserttetidx].quality;
	cavtet.depth = 0;
	cavtet.outfacesize = 0;
	cavtet.parents[0] = NOCAVITYTET;
	cavtet.parents[1] = NOCAVITYTET;
	cavtet.parents[2] = NOCAVITYTET;
	outcavity[outcavitysize++] = cavtet;

	/* initialize candidate face list */
	for (i = 0; i < 4; i++)
	{
		F[Fsize++] = inserttetidx*4+i;
	}

	/* now, as long as we have candidate faces */
	while (Fsize > 0)
	{
		//cuPrintf("\nFaces in F set: \n");
		//for (i = 0; i < Fsize; i++)
		//{
		//	cuPrintf("%d ", F[i]);
		//}
		//cuPrintf("\n");

		if (outcavitysize > MAXCAVITYTETS || Fsize > MAXCAVITYFACES || 
			Wsize > MAXCAVITYFACES || Bsize > MAXCAVITYTETS)
			return;

		/* pull a face out of F */
		currf = F[Fsize-1];
		--Fsize;

		/* get t, the tet on the other side of this face 
		* outin[0]: outward facing tet; outin[1]: inward facing tet*/
		oppf = opposite_halfface(face, halfface, currf);
		if (oppf != NOCAVITYFACE)
			outin[0] = oppf/4;
		else
			outin[0] = NOCAVITYTET;
		outin[1] = currf/4;

		/* the inward facing tet should already be in the output cavity.
		* find it, and add this face as an outgoing face */
		tetidx = findCavityTet(outcavity, outcavitysize, outin[1]);
		if (tetidx == outcavitysize)
			return;

		/* compute the quality of the cavity tet with this face */
		halfFacePoints(point, halfface, currf, fpoint);
		
		cavface.handle = currf;
		cavface.quality = tetquality(pnew, fpoint[0], fpoint[1], fpoint[2], qualmeasure);

		/* check to make sure it's not a ghost tet */
		if (outin[0] == GHOSTTET)
		{
			/* note that this face has no child, and assign it to its parent tet */
			cavface.child = NOCAVITYTET;
			if (outcavity[tetidx].outfacesize < 5)
				outcavity[tetidx].outfaces[(outcavity[tetidx].outfacesize)++] = cavface;

			/* add this face to the wall list */
			W[Wsize++] = currf;
			continue;
		}

		// fetch the outward facing tet
		currtet = outin[0];

		// fetch the other faces of current tet
		int idx = 0;
		j = currtet*4+4;
		for (i = currtet*4; i < j; ++i)
		{
			/* except the opposite half face of the current face*/
			if (i == oppf)
				continue;
			if (idx < 3)
				otherfaces[idx ++] = i;
		}

		/* is t a cavity tet? */
		if (findCavityTet(outcavity, outcavitysize, currtet) != outcavitysize)
		{
			/* we need to add this face to the parent tet, indicating 
			that it has no child because the tet on the other side
			doesn't depend in it's removal to exist */
			cavface.child = NOCAVITYTET;
			if (outcavity[tetidx].outfacesize < 5)
				outcavity[tetidx].outfaces[(outcavity[tetidx].outfacesize)++] = cavface;

			/* yes, do nothing */
			continue;
		}

		/* is t a blocking tet? */
		if (find_array_int(B, Bsize, currtet) != Bsize)
		{
			/* if there is one other wall face of this tet, and the other two faces
			are visible from v, we can add this tet to the cavity. */
			Wcount = 1;
			facing = true;

			j = currtet*4+4;
			for (i = currtet*4; i < j; ++i)
			{
				/* except the current face */
				if (i == oppf)
					continue;

				/* is this face already marked as a wall ? 
				* testing the opposite face which may be visited before */
				if (find_array_int(W, Wsize, opposite_halfface(face, halfface, i)) != Wsize)
				{
					Wcount++;
				}
				else
					/* it's not a wall... is it oriented toward v? */
				{
					halfFacePoints(point, halfface, i, fpoint);
					if (orient(pnew, fpoint[0], fpoint[1], fpoint[2]) <= MINFACING)
					{
						facing = false;
					}
				}
			}

			/* only allow tets with three parents if we are allowing vertex deletion */
			if ((Wcount == 2 || Wcount == 3) && facing)
			{
				/* this tet can be added to the cavity */

				/* remove it from B */
				i = find_array_int(B, Bsize, currtet);
				if (i < Bsize)
				{
					B[i] = B[Bsize-1];
					-- Bsize;
				}
				/* add it to C */
				//C[Csize++] = currtet;

				/* add this tet to the output cavity */
				cavtet.handle = currtet;
				cavtet.outfacesize = 0;
				cavtet.quality = tets[currtet].quality;

				/* we know one parent must be the one we found above */
				cavtet.parents[0] = outcavity[tetidx].handle;
				/* the depth is one more than the parent depth */
				cavtet.depth = outcavity[tetidx].depth + 1;
				/* if this is a new deepest, remember it */
				if (cavtet.depth > deepest) deepest = cavtet.depth;
				outcavity[outcavitysize++] = cavtet;

				/* add this face to the parent tet with the correct child */
				//cavface.child = cavtet;
				if (outcavity[tetidx].outfacesize < 5)
					outcavity[tetidx].outfaces[outcavity[tetidx].outfacesize++] = cavface;

				/* remove any faces that were in W, add others to F. Handle output
				tet faces that need to be added */
				numnonwall = 0;

				/* first, handle all wall face so we can set the correct depth */
				j = currtet*4+4;
				for (i = currtet*4; i<j; ++ i)
				{
					/* except the current face */
					if (i == oppf)
						continue;

					int opphf = opposite_halfface(face, halfface, i);
					/* is this already a wall face? */
					int idx = find_array_int(W, Wsize, opphf);
					if (idx != Wsize)
					{
						W[idx] = W[Wsize-1];
						--Wsize;

						/* because this face was previously a wall face,
						it has some cavity tet that it belongs to. find
						this tet in the output cavity and set it's child face */
						idx = findCavityTet(outcavity, outcavitysize, opphf/4);
						if (idx == outcavitysize)
							return;

						/* add this face to the parent tet's outgoing faces */
						cavface.handle = opphf;
						halfFacePoints(point, halfface, opphf, fpoint);
						cavface.quality = tetquality(pnew, fpoint[0], fpoint[2], fpoint[1], qualmeasure);
						cavface.child = currtet;

						/* make sure that this face is already in this tet */
						if (outcavity[idx].outfacesize < 5)
							outcavity[idx].outfaces[outcavity[idx].outfacesize++] = cavface;

						/* assign the parent tet as the second parent of the new cavity tet */
						if (outcavity[outcavitysize-1].parents[1] == NOCAVITYTET)
						{
							outcavity[outcavitysize-1].parents[1] = outcavity[idx].handle;
						}
						else
						{
							outcavity[outcavitysize-1].parents[2] = outcavity[idx].handle;
						}

						/* if this parent has a lesser depth value, update new tet's depth to be the lesser */
						if (outcavity[idx].depth < outcavity[outcavitysize-1].depth)
						{
							outcavity[outcavitysize-1].depth = outcavity[idx].depth;
						}
					}
					else
					{
						/* record this non-wall face for potential addition to F later */
						nonwall[numnonwall++] = i;
					}
				}

				for (i = 0; i < numnonwall; ++i)
				{
					/* this is a newly-uncovered face. there could be more tets behind it, so
					we should add it to F, if the current tet's depth isn't more than the max */
					if (outcavity[outcavitysize-1].depth < CAVDEPTHLIMIT)
					{
						F[Fsize++] = nonwall[i];
					}
					/* we should artificially make this a wall face so the cavity doesn't get deeper */
					else
					{
						/* construct output face */
						cavface2.handle = nonwall[i];
						halfFacePoints(point, halfface, nonwall[i], fpoint);
						cavface2.quality = tetquality(pnew, fpoint[0], fpoint[2], fpoint[1], qualmeasure);
						cavface2.child = NOCAVITYTET;

						/* add it to parent tet */
						if (outcavity[outcavitysize-1].outfacesize < 5)
							outcavity[outcavitysize-1].outfaces[outcavity[outcavitysize-1].outfacesize++] = cavface2;

						W[Wsize++] = nonwall[i];
					}
				}
			}
			else
			{
				/* note that this face has no child, and assign it to its parent tet */
				cavface.child = NOCAVITYTET;
				if (outcavity[tetidx].outfacesize < 5)
					outcavity[tetidx].outfaces[outcavity[tetidx].outfacesize++] = cavface;

				/* add f to W, it borders a blocking tet */
				W[Wsize++] = currf;
			}
			continue;
		}

		/* t is neither a blocking tet nor a cavity tet */
		/* check to see if the three other faces of the tet are facing v */
		bool check = true;
		for (i = 0; i < 3; ++i)
		{
			/* fetch the face points */
			halfFacePoints(point, halfface, otherfaces[i], fpoint);
			/* the order of points should be reversed */
			if (orient(pnew, fpoint[2], fpoint[1], fpoint[0]) != 1)
			{
				check = false;
				break;
			}
		}
		if (check)
		{
			/* yes! we can add this tet to the cavity */
			//C[Csize++] = currtet;

			/* add this tet to the output cavity */
			cavtet.handle = currtet;
			cavtet.outfacesize = 0;
			cavtet.quality = tets[currtet].quality;

			/* it's parent must be the parent above */
			cavtet.parents[0] = outcavity[tetidx].handle;
			/* depth is one deeper than parent */
			cavtet.depth = outcavity[tetidx].depth + 1;
			/* if this is a new deepest, note it */
			if (cavtet.depth > deepest) deepest = cavtet.depth;
			outcavity[outcavitysize++] = cavtet;

			/* note the current face's child in the parent tet */
			cavface.child = currtet;
			if (outcavity[tetidx].outfacesize < 5)
				outcavity[tetidx].outfaces[outcavity[tetidx].outfacesize++] = cavface;

			/* add t's three (outward oriented) faces to F, if the current tet isn't too deep */
			if (cavtet.depth < CAVDEPTHLIMIT)
			{
				F[Fsize++] = otherfaces[0];
				F[Fsize++] = otherfaces[1];
				F[Fsize++] = otherfaces[2];
			}
			else
			{
				/* construct output face */
				cavface2.child = NOCAVITYTET;
				for (i = 0; i < 3; i++)
				{
					cavface2.handle = otherfaces[i];
					halfFacePoints(point, halfface, otherfaces[i], fpoint);
					cavface2.quality = tetquality(pnew, fpoint[0], fpoint[1], fpoint[2], qualmeasure);
					/* add it to parent tet */
					outcavity[outcavitysize-1].outfaces[outcavity[outcavitysize-1].outfacesize++] = cavface2;
				}
				W[Wsize++] = otherfaces[0];
				W[Wsize++] = otherfaces[1];
				W[Wsize++] = otherfaces[2];
			}
		}
		else
		{
			/* this is a blocking tet, add it to B */
			B[Bsize++] = currtet;

			/* note the current face in the parent tet */
			cavface.child = NOCAVITYTET;
			if (outcavity[tetidx].outfacesize < 5)
			{
				outcavity[tetidx].outfaces[outcavity[tetidx].outfacesize++] = cavface;
			}

			/* add the current face to the wall face list */
			W[Wsize++] = currf;
		}
	}

	/* record the maximum depth */
	//cavdeep = deepest;
}

// insert sort
__device__ void sort_cavityedge(CavityEdge *edges, int edgesize)
{
	int i, j;
	CavityEdge tedge;
	for (i = 1; i < edgesize; ++i)
	{
		tedge = edges[i];
		for (j = i; j > 0; --j)
		{
			if (edges[j-1].qual > tedge.qual)
				edges[j] = edges[j-1];
			else
				break;
		}
		edges[j] = tedge;
	}
}


/* recursively label parents and children as cavity tets */
__device__ void cavityLabel(CavityTet *cavity, int cavitysize, int ctetidx_)
{
	int i;
	int parenttetidx;
	int childtetidx;
	int tetstack[MAXSTACKTET];
	int tetstacksize = 0;
	int ctetidx = ctetidx_;

	/* this tet shouldn't yet be labeled */
	if (cavity[ctetidx].label != NOLABEL)
		return;

	tetstack[tetstacksize++] = ctetidx;

	while (tetstacksize > 0)
	{
		ctetidx = tetstack[tetstacksize-1];
		-- tetstacksize;

		if (cavity[ctetidx].label != NOLABEL)
			continue;
		/* label this tet as in the cavity */

		cavity[ctetidx].label = CAVLABEL;

		/* go through all parents in the original graph */
		for (i = 0; i < 3; ++i)
		{
			if (cavity[ctetidx].parents[i] != NOCAVITYTET)
			{
				/* if this parent is unlabeled, label it */
				parenttetidx = findCavityTet(cavity, cavitysize, cavity[ctetidx].parents[i]);
				if (parenttetidx != cavitysize && cavity[parenttetidx].label == NOLABEL)
				{
					tetstack[tetstacksize++] = parenttetidx;
					//cavityLabel(cavity, cavitysize, parenttetidx);
				}
			}
		}

		/* go through all children in H */
		for (i = 0;  i < cavity[ctetidx].outfacesize; ++i)
		{
			/* check if this edge is in H */
			if (cavity[ctetidx].outfaces[i].inH == true)
			{
				/* this can't be an edge leading to t... we should never add those */
				childtetidx = findCavityTet(cavity, cavitysize, cavity[ctetidx].outfaces[i].child);
				if (childtetidx != cavitysize && cavity[childtetidx].label == NOLABEL)
				{
					tetstack[tetstacksize++] = childtetidx;
					//cavityLabel(cavity, cavitysize, childtetidx);
				}
			}
		}
	}
}

/* recursively label parents and children as anti-cavity tets */
__device__ void antiCavityLabel(CavityTet *cavity, int cavitysize, int ctetidx_)
{
	int parenttetidx;
	int childtetidx;
	int parent,edgetochild;
	int i, j;
	int tetstack[MAXSTACKTET];
	int tetstacksize = 0;
	int ctetidx = ctetidx_;

	/* this tet shouldn't yet be labeled */
	if (cavity[ctetidx].label != NOLABEL)
		return;

	tetstack[tetstacksize++] = ctetidx;

	while (tetstacksize > 0)
	{
		ctetidx = tetstack[tetstacksize-1];
		-- tetstacksize;

		if (cavity[ctetidx].label != NOLABEL)
			continue;
		/* label this tet as in the anticavity */
		cavity[ctetidx].label = ANTICAVLABEL;

		/* go through all parents in H */
		for (i = 0; i < 3; ++i)
		{
			parent = cavity[ctetidx].parents[i];

			if (parent != NOCAVITYTET)
			{
				/* is this parent unlabeled ? */
				parenttetidx = findCavityTet(cavity, cavitysize, cavity[ctetidx].parents[i]);
				if (parenttetidx == cavitysize && cavity[parenttetidx].label != NOLABEL)
				{
					continue;
				}

				/* find the edge from this parent down to the child */
				edgetochild = -1;
				for (j = 0; j < cavity[parenttetidx].outfacesize; j++)
				{
					if (cavity[parenttetidx].outfaces[j].child == cavity[ctetidx].handle)
					{
						edgetochild = j;
					}
				}

				/* is this edge in H? */
				if (cavity[parenttetidx].outfaces[edgetochild].inH == true)
				{
					tetstack[tetstacksize++] = parenttetidx;
					//antiCavityLabel(cavity, cavitysize, parenttetidx);
				}
			}
		}

		/* go through all children in original graph G */
		for (i = 0; i < cavity[ctetidx].outfacesize; i++)
		{
			/* if the child is t, it's the end and is already labeled. move on */
			if (cavity[ctetidx].outfaces[i].child == NOCAVITYTET)
			{
				continue;
			}

			childtetidx = findCavityTet(cavity, cavitysize, cavity[ctetidx].outfaces[i].child);
			if (childtetidx != cavitysize && cavity[childtetidx].label == NOLABEL)
			{
				tetstack[tetstacksize++] = childtetidx;
				//antiCavityLabel(cavity, cavitysize, childtetidx);
			}
		}
	}
}


// smooth the inserting vertex with the cavity faces information
__device__ bool smoothInsertVertex(float *points, cu_halfface *halfface, int *cavityfaces, int cavityfacesize, 
								   float p[3], int qualmetric, float &worstcavity)
{
	float newp[3];
	float qualbefore, qualafter;
	float tempqual;
	int neipoint[MAXCAVITYFACES];
	int neipointcd[MAXCAVITYFACES][3];
	int neipointsize = 0;
	int i, j, k;
	float fpoint[3][3];
	int fpidx;
	//double point[3];
	//double **tetp;
	//int tetcnt;
	//VolumeMesh::Nonsmoother nonsmoother;
	qualbefore = qualafter = 1.0;

	// for each face of cavity
	for (i = 0; i < cavityfacesize; i++)
	{
		// for each point of current face
		for (j = 0; j < 3; j++)
		{
			fpidx = halfface[cavityfaces[i]].pointhandle[j];
			// check if it has been recorded
			for (k = 0; k < neipointsize; ++k)
			{
				if (neipoint[k] == fpidx)
					break;
			}
			// if the point has not been recorded
			if (k == neipointsize)
			{
				neipoint[neipointsize++] = fpidx;
				neipointcd[neipointsize-1][0] = points[fpidx*3];
				neipointcd[neipointsize-1][1] = points[fpidx*3+1];
				neipointcd[neipointsize-1][2] = points[fpidx*3+2];
			}
		}

		halfFacePoints(points, halfface, cavityfaces[i], fpoint);
		tempqual = tetquality(p, fpoint[0], fpoint[1], fpoint[2], qualmetric);
		if (qualbefore > tempqual)
			qualbefore = tempqual;
	}

	/* smooth the point */
	newp[0] = newp[1] = newp[2] = 0.0;
	for (i = 0; i < neipointsize; ++i)
	{
		newp[0] += neipointcd[i][0];
		newp[1] += neipointcd[i][1];
		newp[2] += neipointcd[i][2];
	}
	newp[0] /= neipointsize;
	newp[1] /= neipointsize;
	newp[2] /= neipointsize;
	// fetch smoothing information
	//tetcnt = faces.size()/3;
	//tetp = new double *[tetcnt];
	//for (int i = 0; i < tetcnt; i++)
	//{
	//	tetp[i] = new double[9];
	//	for (int j = 0; j < 3; j++)
	//	{
	//		tetp[i][j*3]   = faces[i*3+j][0];
	//		tetp[i][j*3+1] = faces[i*3+j][1];
	//		tetp[i][j*3+2] = faces[i*3+j][2];
	//	}
	//}
	//point[0] = p[0]; point[1] = p[1]; point[2] = p[2];

	// do vertice smoothing
	//nonsmoother.setTmesh(tmesh);
	//nonsmoother.setQualityKind(MINSINE2);
	//nonsmoother.NonsmoothVertice(point, tetp, tetcnt, point);
	/* end of smoothing */

	// calculate local quality after smoothing
	//newp[0] = point[0];newp[1] = point[1];newp[2] = point[2];
	for (i = 0; i < cavityfacesize; i++)
	{
		halfFacePoints(points, halfface, cavityfaces[i], fpoint);
		tempqual = tetquality(p, fpoint[0], fpoint[1], fpoint[2], qualmetric);
		if (qualafter > tempqual)
			qualafter = tempqual;
	}

	if (qualafter > qualbefore)
	{
		p[0] = newp[0];
		p[1] = newp[1];
		p[2] = newp[2];
		worstcavity = qualafter;
		return true;
	}
	worstcavity = qualbefore;
	return false;
}

/* ALGORITHM 3 + CONSIDER DELETED TETS*/
/* given the DAG representing the largest possible cavity,
find the subcavity that has the lexicographically maximum
quality, then dig it out and connect the exposed faces to 
the new vertex, returning the new tets as well as the worst
quality in this cavity */
__device__ void maxCavity(float *points, cu_halfface *halfface, float pnew[3], CavityTet *cavity, int cavitysize,
						  int *erasetetras, int &erasetetrasize, int *outputfaces, int &outputfacesize, 
						  float &worstdelete, float &cavityqual, int qualmeasure)
{
	//bool foundparent = false;                /* did we find face parent? */
	struct CavityEdge edges[MAXCAVITYEDGES];      /* list of all edges in the DAG */
	int numedges = 0;
	//CavityTet t;                   /* the virtual node t, at the end of array */
	int parentlabel, childlabel;             /* the groups that contain parent and child of an edge */
	float depthfactor;
	//float qual;
	//int deepest = 0;
	int parenttetidx, childtetidx;
	CavityEdge tmpedge;
	int i, j;
	float cavitydepthtable[DEPTHTABLESIZE] = {1.0, 1.6, 2.3, 2.9, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3};

	/* initialize t as a ghost tet */
	//t.handle = NOCAVITYTET;

	/* now, proceed through the DAG, recording all edges, and deleting child
	information to remove the edges, producing D' */
	for (i = 0; i < cavitysize; i++)
	{
		/* initialize this tet's label. set s = cavity, others no label */
		if (i == 0)
		{
			cavity[i].label = CAVLABEL;
		}
		else
		{
			cavity[i].label = NOLABEL;
		}

		/* compute how much to exagerrate quality because of depth */
		if (cavity[i].depth < DEPTHTABLESIZE)
		{
			depthfactor = cavitydepthtable[cavity[i].depth];
		}
		else
		{
			depthfactor = cavitydepthtable[DEPTHTABLESIZE - 1];
		}

		/* for each outgoing face */
		for (j = 0; j < cavity[i].outfacesize; ++j)
		{
			tmpedge.label = EDGELABEL;
			tmpedge.parent = cavity[i].handle;
			tmpedge.qual = cavity[i].outfaces[j].quality * depthfactor;
			/* save which outgoing face this child was */
			tmpedge.childnum = j;

			if (cavity[i].outfaces[j].child == NOCAVITYTET)
			{
				/* if it has no child tet, make up a virtual edge to t */
				tmpedge.child = GHOSTTET;
			}
			else
			{
				/* set edge child to be actual child */
				tmpedge.child = cavity[i].outfaces[j].child;
			}
			edges[numedges++] = tmpedge;

			/* initialize H by setting no edges to be in it */
			cavity[i].outfaces[j].inH = false;
		}
	}

	/* now, sort the edges in order of ascending quality */
	sort_cavityedge(edges, numedges);

	//cuPrintf("\nsorted edges %d :\n", numedges);
	//for (i = 0; i < numedges; i++)
	//{
	//	cuPrintf("edge %d : qual(%f) label(%d) parent(%d) child(%d) childnum(%d)\n", 
	//		     i, edges[i].qual, edges[i].label, edges[i].parent, edges[i].child, edges[i].childnum);
	//}

	/* go through each edge, adding it to D' if it doesn't
	connect s to t */
	for (i = 0; i < numedges; i++)
	{
		/* find parent cavity */
		parenttetidx = findCavityTet(cavity, cavitysize, edges[i].parent);

		/* check parent's label */
		parentlabel = cavity[parenttetidx].label;
		/* check child's label */
		childtetidx = findCavityTet(cavity, cavitysize, edges[i].child);
		if (edges[i].child == GHOSTTET)
			/* this child is labeled automatically as anti-cavity */
			childlabel = ANTICAVLABEL;
		else
			childlabel = cavity[childtetidx].label;

		/* if the parent is in the cavity */
		if (parentlabel == CAVLABEL)
		{
			/* and the child is in the anti-cavity */
			if (childlabel == ANTICAVLABEL)
			{
				/* record output face from parent to child, record the half face of child */
				outputfaces[outputfacesize++] = cavity[parenttetidx].outfaces[edges[i].childnum].handle;
			}
			/* otherwise, if the child isn't labeled */
			else
			{
				if (childlabel == NOLABEL)
				{
					cavityLabel(cavity, cavitysize, childtetidx);
				}
			}
		}
		/* parent isn't labeled cavity. */
		else
		{
			/* is the parent wholly unlabeled ? */
			if (parentlabel == NOLABEL)
			{
				/* is the child labeled anti-cavity ? */
				if (childlabel == ANTICAVLABEL)
				{
					antiCavityLabel(cavity, cavitysize, parenttetidx);
				}
				else
					/* neither the parent nor the child is labeled */
				{
					/* add the edge from parent to child to H */
					cavity[parenttetidx].outfaces[edges[i].childnum].inH = true;
				}
			}
		}
	}

	/* keep track of what the deepest tet in the final cavity was */
	//deepest = 0;

	/* delete all tets labeled as cavity */
	worstdelete = HUGEFLOAT;
	for (i = 0; i < cavitysize; i++)
	{
		if (cavity[i].label == CAVLABEL)
		{
			/* is this the worst quality tet we're deleting? */
			if (cavity[i].quality < worstdelete)
			{
				worstdelete = cavity[i].quality;
			}

			/* is this the deepest tet we've encountered? */
			//if (cavity[i].depth > deepest)
			//{
			//	deepest = cavity[i].depth;
			//}
			erasetetras[erasetetrasize++]= cavity[i].handle;
		}
	}

	// smooth insert point
	smoothInsertVertex(points, halfface, outputfaces, outputfacesize, pnew, qualmeasure, cavityqual);
}



__global__ void vertex_Insertion_Explore(float *points, int pointcnt, cu_InsertTet *tets, int tetcnt, cu_face *face, int facecnt, 
										 cu_halfface *halfface, int halffacecnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	float quality;
	cu_InsertTet currtet;
	float tetpoint[4][3];
	int i;

	while(tid < tetcnt)
	{
		//if (tid != 1653)
		//{
		//	tid += offset;
		//	continue;
		//}
		//if (tid != 1439 && tid != 2263 && tid != 1412)
		//{
		//	tid += offset;
		//	continue;
		//}

		//cuPrintf("\n\ncurrent tid : %d\n", tid);

		currtet = tets[tid];

		// calculate current tet's quality
		int pidx;
		for (i = 0; i < 4; i++)
		{
			pidx = currtet.v[i];
			tetpoint[i][0] = points[3*pidx];
			tetpoint[i][1] = points[3*pidx+1];
			tetpoint[i][2] = points[3*pidx+2];
		}
		quality = tetquality(tetpoint, qualmeasure);
		tets[tid].quality = quality;
		currtet.quality = quality;

		//cuPrintf("\ntet vertex : %d %d %d %d\n", currtet.v[0], currtet.v[1], currtet.v[2], currtet.v[3]);
		//cuPrintf("tet quality : %f\n", currtet.quality);

		// select tetra with bad quality to do vertex insertion
		if (quality > VERTEX_INSERTION_QUALITY_THRESHOLD)
		{
			tid += offset;
			continue;
		}

		//cuPrintf("\nThis tet is a target tet\n");

		int outcavitysize = 0;
		CavityTet outcavity[MAXCAVITYTETS];
		float worstdelete = HUGEFLOAT;
		float cavityqual = HUGEFLOAT;

		/* build the cavity dag */
		float barycenter[3];
		for (i = 0; i < 3; i++)
		{
			barycenter[i] = (tetpoint[0][i] + tetpoint[1][i] + tetpoint[2][i] + tetpoint[3][i])/4.0;
		}

		//cuPrintf("\nCenter point : %f %f %f\n", barycenter[0], barycenter[1], barycenter[2]);

		buildcavitydag(points, tets, face, halfface, barycenter, tid, currtet, tetpoint, qualmeasure, outcavity, outcavitysize);

		//cuPrintf("\nout cavity tet %d :\n", outcavitysize);
		//for (i = 0; i < outcavitysize; i++)
		//{
		//	cuPrintf("%d ", outcavity[i].handle);
		//}
		//cuPrintf("\n");


		/* build the cavity of maximum lexicographic quality */
		int outputfaces[MAXCAVITYFACES];
		int outputfacesize = 0;
		int erasetetras[MAXCAVITYTETS];  /* tetras which are going to be erased */
		int erasetetrasize = 0;
		maxCavity(points, halfface, barycenter, outcavity, outcavitysize, erasetetras, erasetetrasize, 
			      outputfaces, outputfacesize, worstdelete, cavityqual, qualmeasure);

		//cuPrintf("\nerase tet %d :\n", erasetetrasize);
		//for (i = 0; i < erasetetrasize; i++)
		//{
		//	cuPrintf("%d ", erasetetras[i]);
		//}
		//cuPrintf("\n");


		//cuPrintf("\ncavity face %d :\n", outputfacesize);
		//for (i = 0; i < outputfacesize; i++)
		//{
		//	cuPrintf("%d ", outputfaces[i]);
		//}
		//cuPrintf("\n");

		//cuPrintf("\n new point : %f %f %f\n", barycenter[0], barycenter[1], barycenter[2]);

		//cuPrintf("\nquality before: %f     quality after: %f\n", worstdelete, cavityqual);

		/* did we succeed? */
		if (cavityqual > (worstdelete + MINIMPROVEMENT))
		{
			tets[tid].deletetetcnt = erasetetrasize;
			for (i = 0; i < erasetetrasize; i++)
			{
				tets[tid].deletetet[i] = erasetetras[i];
			}
			tets[tid].cavityfacecnt = outputfacesize;
			for (i = 0; i < outputfacesize; i++)
			{
				tets[tid].cavityface[i] = outputfaces[i];
			}
			tets[tid].insertpoint[0] = barycenter[0];
			tets[tid].insertpoint[1] = barycenter[1];
			tets[tid].insertpoint[2] = barycenter[2];
			tets[tid].val = worstdelete;
		}
		tid += offset;
	}
}

__global__ void vertex_insertion_tetquality(float *points, cu_InsertTet *tets, int tetcnt, int qualmeasure)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	float tetpoint[4][3];
	int pidx;

	while(tid < tetcnt)
	{
		// if the tet is invalid
		if (tets[tid].v[0] == -1)
		{
			tid += offset;
			continue;
		}

		// get tetra points
		for (int j = 0; j < 4; j++)
		{
			pidx = tets[tid].v[j];
			tetpoint[j][0] = points[3*pidx];
			tetpoint[j][1] = points[3*pidx+1];
			tetpoint[j][2] = points[3*pidx+2];
		}
		// calculate tetra quality
		tets[tid].quality = tetquality(tetpoint, qualmeasure);

		tid += offset;
	}
}

void inserTet_selection(cu_InsertTet *tets, int tetcnt, int *selecttet, int &selecttetcnt, int &maxincreasetet)
{
	int stetcnt = 0;
	int* stet = new int[tetcnt];
	int i, j;
	int temp;
	maxincreasetet = 0;

	/* pick out the tets by vertex insertion result
	and sort them*/
	for (i = 0; i < tetcnt; i++)
	{
		if (tets[i].deletetetcnt == 0)
			continue;

		// insert the tet
		for (j = stetcnt-1; j > -1 ; j--)
		{
			if (tets[stet[j]].val > tets[i].val)
				stet[j+1] = stet[j];
			else 
				break;
		}
		stet[j+1] = i;
		++ stetcnt;
	}

	/* if there are no any tet meeting the requirement */
	if (!stetcnt)
	{
		selecttetcnt = 0;
		selecttet = NULL;
		return;
	}

	/* select a new tet set in which the cavities of any two of them are not overlapped */
	bool isok;
	bool *tetflag;
	tetflag = new bool[tetcnt];
	memset(tetflag, 1, tetcnt*sizeof(bool));

	selecttet = new int[stetcnt];
	// push the first tet
	selecttet[0] = stet[0];
	selecttetcnt = 1;

	// get increase tet num
	temp = tets[stet[0]].cavityfacecnt - tets[stet[0]].deletetetcnt;
	if (temp > maxincreasetet)
		maxincreasetet = temp;

	// set flags of cavity tets
	for (i = 0; i < tets[stet[0]].deletetetcnt; i++)
	{
		tetflag[tets[stet[0]].deletetet[i]] = 0;
	}

	/* if the cavity tets of current tetra are available, 
	then add the tet into the array*/
	for (i = 1; i < stetcnt; i++)
	{
		isok = true;
		for (j = 0; j < tets[stet[i]].deletetetcnt; j++)
		{
			if (tetflag[tets[stet[i]].deletetet[j]] == 0)
			{
				isok = false;
				break;
			}
		}

		if (isok)
		{
			selecttet[selecttetcnt++] = stet[i];
			for (j = 0; j < tets[stet[i]].deletetetcnt; j++)
			{
				tetflag[tets[stet[i]].deletetet[j]] = 0;
			}

			// get increase tet num
			temp = tets[stet[i]].cavityfacecnt - tets[stet[i]].deletetetcnt;
			if (temp > maxincreasetet)
				maxincreasetet = temp;
		}
	}
}

// 未考虑face和halfface的拓扑修改
__global__ void vertex_insertion(float *points, int pointcnt, cu_InsertTet *tets, cu_halfface *halfface, 
								 int *selecttet, int selecttetcnt, int *increasetets, int maxincreasetet)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	cu_InsertTet currtet;
	int pidx;
	int tetidx;
	int tetnum;
	int i;

	while(tid < selecttetcnt)
	{
		currtet = tets[selecttet[tid]];

		// add new vertex
		pidx = pointcnt + tid;
		points[pidx*3]   = currtet.insertpoint[0];
		points[pidx*3+1] = currtet.insertpoint[1];
		points[pidx*3+2] = currtet.insertpoint[2];

		// add new tet
		if (currtet.deletetetcnt < currtet.cavityfacecnt)
			tetnum = currtet.deletetetcnt;
		else
			tetnum = currtet.cavityfacecnt;

		for (i = 0; i < tetnum; i++)
		{
			tetidx = currtet.deletetet[i];
			tets[tetidx].v[0] = pidx;
			tets[tetidx].v[1] = halfface[currtet.cavityface[i]].pointhandle[0];
			tets[tetidx].v[2] = halfface[currtet.cavityface[i]].pointhandle[1];
			tets[tetidx].v[3] = halfface[currtet.cavityface[i]].pointhandle[2];
		}

		if (currtet.deletetetcnt < currtet.cavityfacecnt)
		{
			for (i = currtet.deletetetcnt; i < currtet.cavityfacecnt; i++)
			{
				tetidx = (maxincreasetet*tid + (i-currtet.deletetetcnt))*4;
				increasetets[tetidx]   = pidx;
				increasetets[tetidx+1] = halfface[currtet.cavityface[i]].pointhandle[0];
				increasetets[tetidx+2] = halfface[currtet.cavityface[i]].pointhandle[1];
				increasetets[tetidx+3] = halfface[currtet.cavityface[i]].pointhandle[2];
			}
		}
		else
		{
			for (i = currtet.cavityfacecnt; i < currtet.deletetetcnt; i++)
			{
				tetidx = currtet.deletetet[i];
				tets[tetidx].v[0] = -1;
				tets[tetidx].v[1] = -1;
				tets[tetidx].v[2] = -1;
				tets[tetidx].v[3] = -1;
			}
		}

		tid += offset;
	}
}

extern "C" void cuda_vertexInsertion(float *&points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
									 int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
									 int &succ, float &time)
{
	//FILE *outfile = NULL;
	//outfile = fopen("F:\\kernel_data_output.txt", "w");

	//fprintf(outfile, "In cuda vertexinsertion\n");

	// create face and halfface structure
	//int *meshtets;
	struct cu_face *cface;
	struct cu_halfface *chalfface;
	struct cu_InsertTet *ctet;
	float qualbefore;
	float qualafter;
	int tetcapacity;
	int facecapacity;
	int halffacecapacity;
	int pointcapacity;

	tetcapacity = int(1.1*tetcnt);
	facecapacity = int(1.1*facecnt);
	halffacecapacity = int(1.1*halffacecnt);
	pointcapacity = int(1.1*pointcnt);

	//meshtets = new int[4*tetcapacity];
	ctet = new struct cu_InsertTet[tetcapacity];
	cface = new struct cu_face[facecapacity];
	chalfface = new struct cu_halfface[halffacecapacity];

	//memcpy(meshtets, meshtets_, 4*tetcnt*sizeof(int));
	for (int i = 0; i < tetcnt; i ++)
	{
		ctet[i].v[0] = meshtets[4*i];
		ctet[i].v[1] = meshtets[4*i+1];
		ctet[i].v[2] = meshtets[4*i+2];
		ctet[i].v[3] = meshtets[4*i+3];
		ctet[i].deletetetcnt = 0;
		memset(ctet[i].deletetet, -1, MAXCAVITYTETS*sizeof(int));
		ctet[i].cavityfacecnt = 0;
		memset(ctet[i].cavityface, -1, MAXCAVITYFACES*sizeof(int));
		ctet[i].quality = 1.0;
	}

	// calculate quality before flipping
	//cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

	// CUDA Parallel
	// 分配设备存储空间
	//int loop = 1;
	float *dev_points;
	struct cu_InsertTet *dev_tets;
	struct cu_face *dev_face/*, *dev_tempface*/;
	struct cu_halfface *dev_halfface/*, *dev_temphalfface*/;
	int *dev_selecttet;
	int *dev_increasetets;

	int *selecttet = NULL;
	int selecttetcnt = 0;
	int *increasetet = NULL;
	int maxincreasetet = 0;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (int i = 0; i < facecnt; i++)
	{
		cface[i].hf[0] = face[2*i];
		cface[i].hf[1] = face[2*i+1];
	}

	for (int i = 0; i < halffacecnt; i++)
	{
		chalfface[i].pointhandle[0] = halfface[4*i];
		chalfface[i].pointhandle[1] = halfface[4*i+1];
		chalfface[i].pointhandle[2] = halfface[4*i+2];
		chalfface[i].face = halfface[4*i+3];
	}


	// 点插入会生成新的四面体、半面和面，因此先预分配一部分空间
	cudaMalloc((void**)&dev_points, 3*pointcapacity*sizeof(float));
	cudaMalloc((void**)&dev_tets, tetcapacity*sizeof(cu_InsertTet));
	cudaMalloc((void**)&dev_face, facecapacity*sizeof(cu_face));
	cudaMalloc((void**)&dev_halfface, halffacecapacity*sizeof(cu_halfface));

	cudaMemcpy(dev_points, points, 3*pointcnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_tets, ctet, tetcnt*sizeof(cu_InsertTet), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_face, cface, facecnt*sizeof(cu_face), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_halfface, chalfface, halffacecnt*sizeof(cu_halfface), cudaMemcpyHostToDevice);

	// vertex insertion explore
	// 调用核函数

	//InitGPUSet();
	//cuPrintInit();

	// step 1 : vertex insertion explore
	//int blocks = imin(tetcnt, BlockPerGrid);
	int blocks = min (tetcnt, BlockPerGrid);
	int threads = min((tetcnt+blocks-1)/blocks, ThreadPerBlock);
	vertex_insertion_tetquality<<<blocks, threads>>>(dev_points, dev_tets, tetcnt, qualmeasure);
	vertex_Insertion_Explore<<<blocks, threads>>>(dev_points, pointcnt, dev_tets, tetcnt, dev_face, facecnt, dev_halfface, halffacecnt, qualmeasure);

	//cudaPrintfDisplay(outfile, false);
	//cudaPrintfEnd(); 
	//fclose(outfile);

	// step 2 : tetra selection
	// get updated tetra data
	cudaMemcpy(ctet, dev_tets, tetcnt*sizeof(cu_InsertTet), cudaMemcpyDeviceToHost);
	inserTet_selection(ctet, tetcnt, selecttet, selecttetcnt, maxincreasetet);

	// step 3 : do parallel vertex insertion
	if (selecttetcnt != 0)
	{
		/** check if point set has enough space
		    if not enlarge the point set*/
		if ((selecttetcnt+pointcnt) > pointcapacity)
		{
			while((selecttetcnt+pointcnt) > pointcapacity)
				pointcapacity = int(1.1*pointcapacity);

			float *dev_temp_points;
			cudaMalloc((void**)&dev_temp_points, 3*pointcnt*sizeof(float));
			cudaMemcpy(dev_temp_points, dev_points, 3*pointcnt*sizeof(float), cudaMemcpyDeviceToDevice);
			cudaFree(dev_points);
			cudaMalloc((void**)&dev_points, 3*pointcapacity*sizeof(float));
			cudaMemcpy(dev_points, dev_temp_points, 3*pointcnt*sizeof(float), cudaMemcpyDeviceToDevice);
			cudaFree(dev_temp_points);
		}

		// do vertex insertion
		cudaMalloc((void**)&dev_selecttet, selecttetcnt*sizeof(int));
		cudaMemcpy(dev_selecttet, selecttet, selecttetcnt*sizeof(int), cudaMemcpyDeviceToDevice);

		// apply space for the increasing tets
		
		cudaMalloc((void**)&dev_increasetets, maxincreasetet*selecttetcnt*4*sizeof(int));
		cudaMemset(dev_increasetets, -1, maxincreasetet*selecttetcnt*4*sizeof(int));

		blocks = min (selecttetcnt, BlockPerGrid);
		threads = min((selecttetcnt+blocks-1)/blocks, ThreadPerBlock);
		vertex_insertion<<<blocks, threads>>>(dev_points, pointcnt, dev_tets, dev_halfface, 
			                                  dev_selecttet, selecttetcnt, dev_increasetets, maxincreasetet);

		pointcnt += selecttetcnt;
		delete [] selecttet;
	}


	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	float elaspsedTime;
	cudaEventElapsedTime(&elaspsedTime, start, stop);

	time = elaspsedTime;

	// calculate quality after flipping
	qualafter = 1.0;
	if (selecttetcnt)
	{
		delete [] points;
		points = new float[pointcnt*3];
		cudaMemcpy(points, dev_points, 3*pointcnt*sizeof(float), cudaMemcpyDeviceToHost);
		
		increasetet = new int [maxincreasetet*selecttetcnt*4];
		cudaMemcpy(increasetet, dev_increasetets, maxincreasetet*selecttetcnt*4*sizeof(int), cudaMemcpyDeviceToHost);
		cuda_tetquality(points, pointcnt, increasetet, maxincreasetet*selecttetcnt, qualmeasure, qualafter);
		delete [] increasetet;
	}
	vertex_insertion_tetquality<<<blocks, (tetcnt+blocks-1)/blocks>>>(dev_points, dev_tets, tetcnt, qualmeasure);
	cudaMemcpy(ctet, dev_tets, tetcnt*sizeof(cu_InsertTet), cudaMemcpyDeviceToHost);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	cudaFree(dev_points);
	cudaFree(dev_tets);
	cudaFree(dev_halfface);
	cudaFree(dev_face);

	for (int i = 0; i < tetcnt; i++)
	{
		if (ctet[i].quality < qualafter)
			qualafter = ctet[i].quality;
	}
	qualbefore = 1.0;
	qualbefore_ = qualbefore;
	qualafter_ = qualafter;
	succ = selecttetcnt;

	delete [] ctet;
	delete [] cface;
	delete [] chalfface;
	//delete [] increasetet;
}
/******************** End of Parallel Vertex Insertion ************************/