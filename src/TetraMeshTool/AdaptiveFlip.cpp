#include <TetraMeshTool/AdaptiveFlip.h>
#include <TetraMeshTool/StarBase.h>
#include <TetraMeshTool/top.h>
namespace VolumeMesh
{
	void AdaptiveFlip::situationObservation(int s[2], std::vector<TetraHandle> &tetvec)
	{
		if (!mesh)
			return;

		Point p[4];
		Point trianglepoint[3], pp;
		int flag;
		int state[3];         // 顶点位置统计 0: outside   1: inside   2: on the boundary
		double qual[3];       // 三种质量值
		TetraMesh::VertexHandle vhandle;             // 三角形情况的顶点
		int situation;        // 0：风筝型 1：三角形  -1: 其他
		double threashold = 0.2;
		TetraMesh::HedronIter h_it;
		TetraMesh::HalfEdgeHandle e1, e2;
		bool issuitable;
		//int s[2];
		s[0] = s[1] = 0;
		for (h_it = mesh->hedrons_begin(); h_it != mesh->hedrons_end(); ++h_it)
		{
			if (mesh->is_boundary(h_it.handle()))
				continue;

			tetPoints(mesh, h_it.handle(), p);
			//qual[0] = qualcalculator.tetquality(p[0], p[1], p[2], p[3], QUAL_RADIUSRATIO);
			qual[1] = qualcalculator.tetquality(p[0], p[1], p[2], p[3], QUAL_MINSINE);
			qual[2] = qualcalculator.tetquality(p[0], p[1], p[2], p[3], QUAL_VLRMS3RATIO);

			situation = 0;
			state[0] = state[1] = state[2] = 0;
			//if (/*qual[0] < threashold && */qual[1] < threashold && qual[2] < threashold)
			{
				/* 确定四面体形状 */
				for (int i = 0; i < 4; i++)
				{
					// 将一点投影到另外三个点所在的平面
					pp = p[i];
					for (int j = 0; j < 3; j++)
						trianglepoint[j] = p[(i+j+1)%4];

					vertexTriangleDistance(trianglepoint, pp, flag);
					if (flag == 0)
						++state[0];
					else if (flag == 1)
					{
						++state[1];
						vhandle = TetraMesh::VertexHandle((h_it.handle()<<2)+i);    // 顶点编号记录
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
					TetraMesh::HalfFaceHandle hf_;
					TetraMesh::FaceHalfedgeIter fe_it;
					TetraMesh::HalfEdgeHandle he_[6];
					double edgelength[6];
					int eidx;

					eidx = 0;
					hf_ = TetraMesh::HalfFaceHandle(h_it.handle().idx()<<2);  //四面体的第一个半面
					for (fe_it = mesh->face_half_edge_iter(hf_); fe_it && eidx < 3; ++fe_it)
					{
						// current halfedge
						he_[eidx] = fe_it.handle();
						// fetch opposite halfedge of the current halfedge
						he_[eidx+3] = mesh->prev_half_edge_handle(mesh->mate_half_edge_handle(mesh->next_half_edge_handle(he_[eidx])));
						++eidx;
					}

					for (int eidx = 0; eidx < 3; eidx++)
					{
						e1 = he_[eidx];
						e2 = he_[(eidx+3)%6];

						issuitable = kitesituationcheck(h_it.handle(), e1, e2);
						if (issuitable)
						{
							kiteAdaptiveFlip(h_it.handle(), e1, e2);
							tetvec.push_back(h_it.handle());
							++s[0];
							break;
						}
					}
				}
				else if(situation == 1)//三角型
				{
					issuitable = trianglesituationcheck(h_it.handle(), vhandle);
					if (issuitable)
					{
						triangleAdaptiveFlip(h_it.handle(), vhandle);
						tetvec.push_back(h_it.handle());
						++s[1];
					}
				}
			}
		}
	}

	/* 确认三角型四面体周围的四面体分布是否复合要求 */
	bool AdaptiveFlip::trianglesituationcheck(TetraHandle th_, TetraMesh::VertexHandle vh_)
	{
		// 与vh_周围三个面相邻的三个四面体相交于一点
		TetraMesh::HalfFaceHandle hf[3];
		hf[0] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+1)%4);
		hf[1] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+2)%4);
		hf[2] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+3)%4);

		TetraMesh::HalfFaceHandle ophf[3];
		TetraMesh::PointHandle ph[3];
		for (int i = 0; i < 3; i ++)
		{
			ophf[i] = mesh->opposite_half_face_handle(hf[i]);
			ph[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf[i].idx())); 
		}

		if (ph[0] == ph[1] && ph[0] == ph[2])
			return true;
		return false;
	}

	/* 确认风筝型四面体周围的四面体分布是否复合要求 */
	bool AdaptiveFlip::kitesituationcheck(TetraHandle th_, TetraMesh::HalfEdgeHandle e1, TetraMesh::HalfEdgeHandle e2)
	{
		TetraMesh::HalfFaceHandle hf1[2];
		TetraMesh::HalfFaceHandle hf2[2];
		TetraMesh::HalfFaceHandle ophf1[2];
		TetraMesh::HalfFaceHandle ophf2[2];
		PointHandle p1[2], p2[2];

		hf1[0] = mesh->handle_to_entity(e1).half_face_handle();
		hf1[1] = mesh->handle_to_entity(mesh->mate_half_edge_handle(e1)).half_face_handle();

		hf2[0] = mesh->handle_to_entity(e2).half_face_handle();
		hf2[1] = mesh->handle_to_entity(mesh->mate_half_edge_handle(e2)).half_face_handle();

		for (int i = 0; i < 2; i ++)
		{
			ophf1[i] = mesh->opposite_half_face_handle(hf1[i]);
			p1[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf1[i].idx()));
			ophf2[i] = mesh->opposite_half_face_handle(hf2[i]);
			p2[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf2[i].idx()));
		}
		
		if (p1[0] == p1[1] && p2[0] == p2[1])
			return true;
		return false;
	}

	/* 查找风筝型四面体中的相交的对边 */
	void AdaptiveFlip::findcrossedges(TetraMesh::HedronHandle th_, TetraMesh::HalfEdgeHandle &e1, TetraMesh::HalfEdgeHandle &e2)
	{
		TetraMesh::HalfFaceHandle hf_;
		TetraMesh::FaceHalfedgeIter fe_it;
		TetraMesh::HalfEdgeHandle he_[6];
		double edgelength[6];
		int eidx;

		hf_ = TetraMesh::HalfFaceHandle(th_.idx()<<2);  //四面体的第一个半面

		eidx = 0;
		for (fe_it = mesh->face_half_edge_iter(hf_); fe_it && eidx < 3; ++fe_it)
		{
			// current halfedge
			he_[eidx] = fe_it.handle();
			edgelength[eidx] = (mesh->point(mesh->from_vertex_handle(he_[eidx])) - mesh->point(mesh->to_vertex_handle(he_[eidx]))).norm();
			// fetch opposite halfedge of the current halfedge
			he_[eidx+3] = mesh->prev_half_edge_handle(mesh->mate_half_edge_handle(mesh->next_half_edge_handle(he_[eidx])));
			edgelength[eidx+3] = (mesh->point(mesh->from_vertex_handle(he_[eidx+3])) - mesh->point(mesh->to_vertex_handle(he_[eidx+3]))).norm();
			++eidx;
		}

		// 相交对边中，一定有一条边是四面体中最长的边
		eidx = 0;
		for (int i = 0; i < 6; i ++)
		{
			if (edgelength[i] > edgelength[eidx])
				eidx = i;
		}

		e1 = he_[eidx];
		e2 = he_[(eidx+3)%6];
	}


	/** 
	* the distance between the triangle and a point
	**/
	double AdaptiveFlip::vertexTriangleDistance(Point triangle[3], Point p, int & pStatus)
	{
		double dis;
		dis = 0;

		TetraMesh::Normal edge1 = triangle[1] - triangle[0];
		TetraMesh::Normal edge2 = triangle[2] - triangle[0];
		TetraMesh::Normal faceNormal;
		faceNormal = (edge1 % edge2).normalize();
		pStatus = 0;

		dis = (faceNormal[0] * p[0] + faceNormal[1] * p[1] +
			faceNormal[2] * p[2] - faceNormal[0] * triangle[0][0] -
			faceNormal[1] * triangle[0][1] - faceNormal[2] * triangle[0][2]) /
			faceNormal.norm();

		Point pp;
		pp = p - faceNormal * dis; // the projective point

		pStatus = checkPointInsideTriangle(triangle, pp);
		return dis;
	}

	/** 
	* check if the point is in the triangle
	* return: 1 the point is in the triangle
	          0 the point is out of the triangle
			 -1 the point is on the boundary of the triangle
	**/
	int AdaptiveFlip::checkPointInsideTriangle(Point triangle[3], Point p)
	{
		TetraMesh::Normal edge1 = triangle[1] - triangle[0];
		TetraMesh::Normal edge2 = triangle[2] - triangle[0];
		TetraMesh::Normal faceNormal;
		faceNormal = (edge1 % edge2).normalize();
		TetraMesh::Normal v[3];
		v[0] = (triangle[1] - triangle[0]) % (p - triangle[0]);
		v[1] = (triangle[2] - triangle[1]) % (p - triangle[1]);
		v[2] = (triangle[0] - triangle[2]) % (p - triangle[2]);
		double cosAngle;

		for (int i = 0; i < 3; i ++)
		{
			if (v[i].norm() > EPS)
			{
				cosAngle = (v[i] | faceNormal) / (double)v[i].norm();
			}
			else
			{
				cosAngle = 0;
			}
			if (cosAngle < - EPS)
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


	void AdaptiveFlip::kiteAdaptiveFlip(TetraHandle th_, TetraMesh::HalfEdgeHandle e1, TetraMesh::HalfEdgeHandle e2)
	{
		if (!mesh)
			return;

		if (mesh->is_boundary(th_))
			return;

		TetraMesh::HalfFaceHandle hf1[2];
		TetraMesh::HalfFaceHandle hf2[2];
		TetraMesh::HalfFaceHandle ophf1[2];
		TetraMesh::HalfFaceHandle ophf2[2];
		TetraMesh::HedronHandle tet1[2];
		TetraMesh::HedronHandle tet2[2];
		PointHandle p1[2], p2[2];
		PointHandle endp1[2], endp2[2];
		Point newtet[16];
		double qualtet1[2], qualtet2[2], tetqual;
		double qualnewtet[2];
		double minqualbefore, minqualafter1, minqualafter2, qual;
		Point tetp[4];

		hf1[0] = mesh->handle_to_entity(e1).half_face_handle();
		hf1[1] = mesh->handle_to_entity(mesh->mate_half_edge_handle(e1)).half_face_handle();

		hf2[0] = mesh->handle_to_entity(e2).half_face_handle();
		hf2[1] = mesh->handle_to_entity(mesh->mate_half_edge_handle(e2)).half_face_handle();

		endp1[0] = mesh->point_handle(mesh->from_vertex_handle(e1));
		endp1[1] = mesh->point_handle(mesh->to_vertex_handle(e1));

		endp2[0] = mesh->point_handle(mesh->from_vertex_handle(e2));
		endp2[1] = mesh->point_handle(mesh->to_vertex_handle(e2));

		for (int i = 0; i < 2; i ++)
		{
			ophf1[i] = mesh->opposite_half_face_handle(hf1[i]);
			tet1[i] = mesh->handle_to_entity(ophf1[i]).hedron_handle();
			p1[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf1[i].idx()));
			ophf2[i] = mesh->opposite_half_face_handle(hf2[i]);
			tet2[i] = mesh->handle_to_entity(ophf2[i]).hedron_handle();
			p2[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf2[i].idx()));
		}

		if (p1[0] != p1[1] || p2[0] != p2[1])
			return;

		// 计算原来的质量值
		tetPoints(mesh, th_, tetp);
		tetqual = qualcalculator.tetquality(tetp[0], tetp[1], tetp[2], tetp[3], QUAL_MINSINE);
		minqualbefore = tetqual;
		for (int i = 0; i < 2; i ++)
		{
			tetPoints(mesh, tet1[i], tetp);
			qualtet1[i] = qualcalculator.tetquality(tetp[0], tetp[1], tetp[2], tetp[3], QUAL_MINSINE);
			tetPoints(mesh, tet2[i], tetp);
			qualtet2[i] = qualcalculator.tetquality(tetp[0], tetp[1], tetp[2], tetp[3], QUAL_MINSINE);
			if (minqualbefore > qualtet1[i])
				minqualbefore = qualtet1[i];
			if (minqualbefore > qualtet2[i])
				minqualbefore = qualtet2[i];
		}

		// 组合策略一：e1周围两个四面体+新的两个四面体
		newtet[0] = mesh->point(p2[0]);
		newtet[1] = mesh->point(endp2[0]);
		newtet[2] = mesh->point(endp1[1]);
		newtet[3] = mesh->point(endp1[0]);

		newtet[4] = mesh->point(p2[0]);
		newtet[5] = mesh->point(endp2[1]);
		newtet[6] = mesh->point(endp1[0]);
		newtet[7] = mesh->point(endp1[1]);

		for (int i = 0; i < 2; i ++)
		{
			qualnewtet[i] = qualcalculator.tetquality(newtet[i*4], newtet[i*4+1], 
				                                      newtet[i*4+2], newtet[i*4+3], QUAL_MINSINE);
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
		newtet[8] = mesh->point(p1[0]);
		newtet[9] = mesh->point(endp2[1]);
		newtet[10] = mesh->point(endp2[0]);
		newtet[11] = mesh->point(endp1[0]);

		newtet[12] = mesh->point(p1[0]);
		newtet[13] = mesh->point(endp2[0]);
		newtet[14] = mesh->point(endp2[1]);
		newtet[15] = mesh->point(endp1[1]);

		for (int i = 0; i < 2; i ++)
		{
			qualnewtet[i] = qualcalculator.tetquality(newtet[(i+2)*4], newtet[(i+2)*4+1], 
				                                      newtet[(i+2)*4+2], newtet[(i+2)*4+3], QUAL_MINSINE);
		}

		minqualafter2 = 1.0;
		for (int i = 0; i < 2; i ++)
		{
			if (minqualafter2 > qualtet2[i])
				minqualafter2 = qualtet2[i];
			if (minqualafter2 > qualnewtet[i])
				minqualafter2 = qualnewtet[i];
		}

		if (minqualafter2 > minqualafter1 && minqualafter2 > minqualbefore)
		{
			mesh->erase_tetrahedron(th_);
			mesh->erase_tetrahedron(tet1[0]);
			mesh->erase_tetrahedron(tet1[1]);
			for (int i = 0; i < 2; i ++)
			{
				mesh->add_tetrahedron(newtet[(i+2)*4], newtet[(i+2)*4+1], 
					                  newtet[(i+2)*4+2], newtet[(i+2)*4+3], &tet1[i]);
			}
		}
		else if (minqualafter1 > minqualafter2 && minqualafter1 > minqualbefore)
		{
			mesh->erase_tetrahedron(th_);
			mesh->erase_tetrahedron(tet2[0]);
			mesh->erase_tetrahedron(tet2[1]);
			for (int i = 0; i < 2; i ++)
			{
				mesh->add_tetrahedron(newtet[i*4], newtet[i*4+1], 
					                  newtet[i*4+2], newtet[i*4+3], &tet2[i]);
			}
		}
	}

	void AdaptiveFlip::triangleAdaptiveFlip(TetraHandle th_, TetraMesh::VertexHandle vh_)
	{
		// 与vh_周围三个面相邻的三个四面体相交于一点
		TetraMesh::HalfFaceHandle hf[3];
		TetraMesh::HalfFaceHandle hf_;
		hf[0] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+1)%4);
		hf[1] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+2)%4);
		hf[2] = TetraMesh::HalfFaceHandle((th_<<2)+(vh_.idx()+3)%4);
		hf_ = TetraMesh::HalfFaceHandle(vh_.idx());

		TetraMesh::HalfFaceHandle ophf[3];
		TetraMesh::PointHandle ph[3];
		TetraMesh::HedronHandle tet[3];
		for (int i = 0; i < 3; i ++)
		{
			ophf[i] = mesh->opposite_half_face_handle(hf[i]);
			tet[i] = mesh->handle_to_entity(ophf[i]).hedron_handle();
			ph[i] = mesh->point_handle(TetraMesh::VertexHandle(ophf[i].idx())); 
		}

		if (ph[0] != ph[1] || ph[0] != ph[2])
			return;

		double qual, minqualbefore, minqualafter;
		std::vector<TetraHandle> tetvec;
		tetvec.push_back(th_);
		for (int i  = 0; i < 3; i ++)
			tetvec.push_back(tet[i]);
		minqualbefore = qualcalculator.minstackquality(mesh, tetvec, QUAL_MINSINE);

		Point fp[3];
		TetraMesh::HalfFaceVertexIter fv_it;
		fv_it = mesh->half_face_vertex_iter(hf_);
		fp[0] = mesh->point(fv_it.handle()); ++fv_it;
		fp[1] = mesh->point(fv_it.handle()); ++fv_it;
		fp[2] = mesh->point(fv_it.handle());

		minqualafter = qualcalculator.tetquality(mesh->point(ph[0]), fp[0], fp[1], fp[2], QUAL_MINSINE);

		if (minqualafter > minqualbefore)
		{
			mesh->erase_tetrahedron(th_);
			for (int i = 0; i < 3; i ++)
			{
				mesh->erase_tetrahedron(tet[i]);
			}
			mesh->add_tetrahedron(mesh->point(ph[0]), fp[0], fp[1], fp[2], &tet[0]);
		}
	}
}