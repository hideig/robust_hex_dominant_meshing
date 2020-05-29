#include <TetraMeshTool/CTetraMeshImprove.h>
#include <fstream>
#include <time.h>
#include <Windows.h>

#define HUGEQUAL 1.0e38
#define MINFLIPIMPROVE 1.0e-5
#define imin(a,b) (a<b?a:b)
namespace VolumeMesh
{
	//void cuda_vertexSmoothing(float *points, int pointcnt, int *hneighbour, int *neighbourcnt, int *pointgroup, int pointgroupcnt, 
	//	float *newpoints, int *hincidenttet, int *incidenttetcnt, int *meshtets, int tetcnt, 
	//	int largestn, int largesttet, int smoothpasscnt, float &time)
	//{

	//}

	//void cuda_tetquality(float *points, int pointcnt, int *meshtets, int tetcnt, int qualmeasure, float &minqual)
	//{

	//}

	//void cuda_flip23(float *points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
	//	int *halfface, int halffacecnt, int qualmeasure, float& qualbefore, float& qualafter, 
	//	int &flipsucc, float &time)
	//{

	//}

	//void cuda_flip32(float *points, int pointcnt, int *meshtets_, int tetcnt, int *flipedge, int edgecnt, 
	//	int qualmeasure, float& qualbefore_, float& qualafter_, int &flipsucc, float &time)
	//{

	//}

	//void cuda_flip32_new(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
	//	int *halfedge, int halfedgecnt, int qualmeasure, float &qualbefore_, float &qualafeter_,
	//	int &flipsucc, float &time)
	//{

	//}

	//void cuda_newflip(float *points, int pointcnt, int *meshtets, int tetcnt, int *face, int facecnt, 
	//	int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
	//	int &flipsucc, float &time)
	//{

	//}

	//void cuda_edgeContraction(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
	//	int *halfedge, int halfedgecnt, int *incidenttet, int* incidenttetcnt, int &largesttetcnt, 
	//	int qualmeasure, float &qualbefore_, float &qualafeter_, int &succ, float &time)
	//{

	//}

	//void cuda_vertexInsertion(float *&points, int pointcnt, int *meshtets_, int tetcnt, int *face, int facecnt, 
	//	int *halfface, int halffacecnt, int qualmeasure, float& qualbefore_, float& qualafter_, 
	//	int &succ, float &time)
	//{

	//}


	void CTetraMeshImprove::pointsandneighbours(float* points, int pointcnt, int** &neighbours, int* neighbourcnt, 
		                                        int** &incidenttet, int* incidenttetcnt, float &time)
	{
		if (mesh == NULL)
			return;

		LARGE_INTEGER li_start, li_finish;
		std::vector<TetraHandle> pointstar;
		std::set<int> pointneibour;
		TetraMesh::PointIter p_it;
		TetraMesh::HedronVertexIter hv_it;
		Point tp;
		PointHandle tph;

		QueryPerformanceCounter(&li_start);
		int pidx = 0;
		// for every point : fetch its coordinate and its neighbour points
		for (p_it = mesh->points_begin(); p_it != mesh->points_end(); ++p_it)
		{
			// get points coordinates
			tp = mesh->point(p_it);
			points[pidx*3]   = (float)tp[0];
			points[pidx*3+1] = (float)tp[1];
			points[pidx*3+2] = (float)tp[2];

			pointneibour.clear();
			pointstar.clear();
			mesh->point_star(p_it.handle(), pointstar);

			incidenttet[pidx] = new int [pointstar.size()];
			incidenttetcnt[pidx] = pointstar.size();
			// get incident tetras & points neighbours
			for (unsigned int i = 0; i < pointstar.size(); i++)
			{
				// get incident tetra
				incidenttet[pidx][i] = pointstar[i].idx();
				// get points neighbours
				for (hv_it = mesh->hedron_vertex_iter(pointstar[i]); hv_it; ++hv_it)
				{
					tph = mesh->point_handle(hv_it.handle());
					if (tph != p_it.handle())
					{
						pointneibour.insert(tph.idx());
					}
				}
			}
			neighbours[pidx] = new int[pointneibour.size()];
			neighbourcnt[pidx] = pointneibour.size();
			copy(pointneibour.begin(), pointneibour.end(), neighbours[pidx]);

			++pidx;
		}
		QueryPerformanceCounter(&li_finish);
		time = 0.0;
		time += float(li_finish.QuadPart - li_start.QuadPart);
	}


	void CTetraMeshImprove::pointsNeighbourTetra(TetraMesh *mesh, int** &incidenttet_, int* &incidenttetcnt, int &largesttetcnt, float &time)
	{
		if (mesh == NULL)
			return;

		LARGE_INTEGER m_liPerfFreq={0}; 
		LARGE_INTEGER li_start, li_finish;
		std::vector<TetraHandle> pointstar;
		TetraMesh::PointIter p_it;
		int **incidenttet;
		int pointcnt;
		int pointstarcnt;
		pointcnt = mesh->size_point();
		incidenttetcnt = new int[pointcnt];
		incidenttet = new int *[pointcnt];
		largesttetcnt = 0;

		QueryPerformanceFrequency(&m_liPerfFreq);  //获取每秒多少CPU Performance Tick  
		QueryPerformanceCounter(&li_start);
		int pidx = 0;
		// for every point : fetch its coordinate and its neighbour points
		for (p_it = mesh->points_begin(); p_it != mesh->points_end(); ++p_it)
		{
			pointstar.clear();
			mesh->point_star(p_it.handle(), pointstar);

			pointstarcnt = pointstar.size();
			incidenttet[pidx] = new int [pointstarcnt];
			incidenttetcnt[pidx] = pointstarcnt;

			if (largesttetcnt < pointstarcnt)
				largesttetcnt = pointstarcnt;

			// get incident tetras
			for (unsigned int i = 0; i < pointstar.size(); i++)
			{
				// get incident tetra
				incidenttet[pidx][i] = pointstar[i].idx();
			}
			++pidx;
		}

		for (int i = 0; i < pointcnt; i ++)
		{
			incidenttet_[i] = new int[largesttetcnt];
			for (int j= 0; j < incidenttetcnt[i]; j++)
			{
				incidenttet_[i][j] = incidenttet[i][j];
			}
		}

		QueryPerformanceCounter(&li_finish);

		time = 0.0;
		time += float(li_finish.QuadPart - li_start.QuadPart);
		time = time * 1000 / m_liPerfFreq.QuadPart;
	}


	void CTetraMeshImprove::getcudapoints(TetraMesh *mesh, float* points)
	{
		if (mesh == NULL && points != NULL)
			return;

		TetraMesh::PointIter p_it;
		Point tp;

		int pidx = 0;
		// for every point : fetch its coordinate and its neighbour points
		for (p_it = mesh->points_begin(); p_it != mesh->points_end(); ++p_it)
		{
			// get points coordinates
			tp = mesh->point(p_it);
			points[pidx*3]   = (float)tp[0];
			points[pidx*3+1] = (float)tp[1];
			points[pidx*3+2] = (float)tp[2];
			++pidx;
		}
	}

	/* quality calculate */
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

	void vector_cross(float *a, float *b, float *c)
	{
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}

	float vector_dot(float *a, float *b)
	{
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

	void vector_add(float *a, float *b, float *c)
	{
		c[0] = a[0] + b[0];
		c[1] = a[1] + b[1];
		c[2] = a[2] + b[2];
	}

	void vector_minus(float *a, float *b, float *c)
	{
		c[0] = a[0] - b[0];
		c[1] = a[1] - b[1];
		c[2] = a[2] - b[2];
	}

	float minsine(float point[4][3])
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

	float tetquality(float point[4][3], int qualmeasure)
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

	float mintetquality(float *points, int *meshtets, int tetcnt, int qualmeasure)
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

	float minstackquality(float *points, int *meshtets, int *tethandlestack, int stackcnt, int qualmeasure)
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

	float minstackquality(float *tetrapoints, int tetracnt, int qualmeasure)
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
	/* end of quality calculate */

	// testing
	void edgeSelecting_new_test(cu_edge *cuedge, int edgecnt, int tetcnt, int *selectedge, int &selectedgecnt)
	{
		int sedgecnt = 0;
		int* sedge = new int[edgecnt];
		int i, j;
		float threashold = 1.0e-5;

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

	// not tested yet
	void getECNewTetras_test(int *newtetra, int *newtetracnt, int *fromtetra, int fromtetracnt, int *totetra, int totetracnt,
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

	// not tested yet
	void edge_contraction_explore_test(float *points, int *meshtets, cu_edge *edge, int edgecnt, cu_halfedge *halfedge, 
		int *incidenttet, int *incidenttetcnt, int largesttetcnt, int qualmeasure)
	{
		int tid = 0;
		int offset = 1;
		int fromv, tov, tempv;
		float fp[3], tp[3], mp[3];
		int edgestar[50];
		int newtetra[200];
		int fromtetra[100];
		int totetra[100];
		int edgestarcnt;
		int newtetracnt;
		int fromtetracnt;
		int totetracnt;
		float newtetrapoint[2400];
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
			edgestarcnt = curredge.halfedgecnt;
			for (i = 0; i < edgestarcnt; i++)
				edgestar[i] = (phe[i])/12;

			// get new tetras
			fromv = halfedge[phe[0]].fromv;
			tov = halfedge[phe[0]].tov;
			fromtetracnt = incidenttetcnt[fromv];
			totetracnt = incidenttetcnt[tov];

			for (i = 0; i < fromtetracnt; i++)
				fromtetra[i] = incidenttet[fromv*largesttetcnt + i];
			for (i = 0; i < totetracnt; i++)
				totetra[i] = incidenttet[tov*largesttetcnt + i];

			getECNewTetras_test(newtetra, &newtetracnt, fromtetra, fromtetracnt, totetra, totetracnt, edgestar, edgestarcnt);

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

	void edge_contraction_test(float *points, int *meshtets, cu_edge *edge, cu_halfedge *halfedge,
		int *selectedge, int selectedgecnt)
	{
		int tid = 0;
		int offset = 1;

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

	void CTetraMeshImprove::cuda_edgeContraction_test(float *points, int pointcnt, int *meshtets, int tetcnt, int *edge, int edgecnt,
		int *halfedge, int halfedgecnt, int *incidenttet, int* incidenttetcnt, int &largesttetcnt, 
		int qualmeasure, float &qualbefore_, float &qualafter_, int &succ_, float &time)
	{
		FILE *outfile = NULL;
		outfile = fopen("F:\\kernel_data_output.txt", "w");

		fprintf(outfile, "My application %d \n", 5);

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
			cedge[i].halfedge = new int[2*cedge[i].halfedgecnt];
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

		float qualbefore, qualafter;
		int *selectedge;
		int selectedgecnt;
		int loop = 1;

		// calculate quality before edge contraction
		cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualbefore);

		int succ = 0;
		while(loop)
		{

			/***********do parallel edge contraction**********/
			// edge contraction explore
			edge_contraction_explore_test(points, meshtets, cedge, edgecnt, chalfedge, incidenttet, incidenttetcnt, largesttetcnt, qualmeasure);

			// edge selecting
			selectedge = new int[edgecnt];
			edgeSelecting_new_test(cedge, edgecnt, tetcnt, selectedge, selectedgecnt);

			// edge contraction (set -1 to the meshtets[] of invalid tets)
			if (selectedgecnt)
			{
				// do edge contraction
				edge_contraction_test(points, meshtets, cedge, chalfedge, selectedge, selectedgecnt);
				succ += selectedgecnt;
			}
			-- loop;
		}

		// calculate quality after flipping
		qualafter = 1.0;
		// calculate quality before flipping
		cuda_tetquality(points, pointcnt, meshtets, tetcnt, qualmeasure, qualafter);

		qualbefore_ = qualbefore;
		qualafter_ = qualafter;
		succ_ = succ;

		delete [] cedge;
		delete [] chalfedge;


	}
}