//#include <TetraMeshTool/top.h>
//#include <stack>
//#include <TetraMeshTool/StarBase.h>
//#include <iterator>
//#include <set>
#include <TetraMeshTool/VertexNonsmooth.h>
#include <fstream>
#include <iostream>
#include <set>
//


namespace VolumeMesh
{
	// 从mesh中读入点，并为各自的数据结构赋值，形成原始拓扑。
	void Nonsmoother::initialinformation()
	{
		TetraMesh::HedronIter h_iter;
		TetraMesh::HedronVertexIter hv_iter;
		TetraMesh::VertexHedronIter vh_iter;
		TetraMesh::PointHedronIter ph_iter;
		TetraMesh::PointIter p_iter;
		TetraMesh::PointHandle ph, ph_2;
		TetraMesh::Point point;
		TetraMesh::VertexHandle vh_1, vh_2;
		int index, point_index;
		int i;
		int incident_hedron_num;
		vertex_num_ = tmesh_->n_point();
		hedron_num_ = tmesh_->size_tetrahedron();
		points_ = new NonsmoothPoint[vertex_num_];//建立对应大小的点数组
		old_points_ = new NonsmoothPoint[vertex_num_];//用于保存最原始的点坐标，用于smooth后的比较使用
		tetrahedrons2_ = new OneTetrahedron2[hedron_num_];//建立对应大小的四面体的数组，每个四面体只存储四个顶点在points_中的index
		incident_tetradrons2_ = new IncidentTetrahedrons2[vertex_num_];//记录下每个点相邻的四面体的index
		for (p_iter = tmesh_->points_begin(); p_iter != tmesh_->points_end(); ++ p_iter)
		{
			ph = p_iter.handle();
			index = ph.idx();
			point = tmesh_->point(ph);
			points_[index].point_coord[0] = point[0];
			points_[index].point_coord[1] = point[1];
			points_[index].point_coord[2] = point[2];
			old_points_[index].point_coord[0] = point[0];
			old_points_[index].point_coord[1] = point[1];
			old_points_[index].point_coord[2] = point[2];
			//记录边界情况
			if(tmesh_->is_boundary(ph))
			{
				points_[index].boundary = true;
				old_points_[index].boundary = true;
			}
			else
			{
				points_[index].boundary = false;
				old_points_[index].boundary = false;
			}
		}
		//clear incident tetrahedron number for every vertices
		for(int i = 0; i < vertex_num_; ++ i)
		{
			incident_tetradrons2_[i].incident_hedron_num = 0;
		}
		//record the tetrahedron
		for (h_iter = tmesh_->hedrons_begin(); h_iter != tmesh_->hedrons_end(); ++ h_iter)
		{
			TetraMesh::HedronHandle hedron;
			OneTetrahedron  temp_hedron;
			hedron = h_iter.handle();
			index = hedron.idx();
			i = 0;
			//设置四面体边界与否信息
			if(tmesh_->is_boundary(hedron))
				tetrahedrons2_[index].boundary = true;
			else
				tetrahedrons2_[index].boundary = false;
			//遍历每一个四面体，将四个顶点的信息存储下来，并对四个顶点所在的incident_tetradron2_中加入该四面体
			for (hv_iter = tmesh_->hedron_vertex_iter(hedron); hv_iter; ++ hv_iter, ++ i)
			{
				vh_1 = hv_iter.handle();
				ph = tmesh_->point_handle(vh_1);
				point_index = ph.idx();
				tetrahedrons2_[index].four_points_index[i] = point_index;
				incident_hedron_num = incident_tetradrons2_[point_index].incident_hedron_num;
				incident_tetradrons2_[point_index].incident_tetrahedrons[incident_hedron_num] = index;
				++ incident_tetradrons2_[point_index].incident_hedron_num;
			}
		}	
	}
	//获得整个model的最小质量，只要遍历每一个四面体分别计算制定的质量即可。
	double Nonsmoother::getMinQuality()
	{
		NonsmoothPoint p[4];
		int point_index;
		double minquality = 1.0;
		double newquality;
		int i, j;
		for (i = 0; i < hedron_num_; ++ i)
		{
			//if (!tetrahedrons2_[i].boundary)
			{
				for (j = 0; j < 4; ++ j)
				{
					point_index = tetrahedrons2_[i].four_points_index[j];
					p[j] = points_[point_index];
				}
				newquality = nonsmoothquality_.getQuality(p, quality_kinds_);
				if(newquality < minquality)
					minquality = newquality;
			}	
		}
		return minquality;
	}
	//获得面的法向，面是p2,p3,p4组成的，其中p1,p2,p3,p4在同一个四面体内部，用叉乘求原始法向，然后根据p1来调整
	//法向的方向，此处有可能因为精度而出现问题
	NonsmoothPoint Nonsmoother::getFacenormal(NonsmoothPoint p1, NonsmoothPoint p2, NonsmoothPoint p3, NonsmoothPoint p4)
	{
		NonsmoothPoint face_normal;
		face_normal = (p2 - p4) % (p3 - p4);
		if(((p1 - p4) | face_normal) >= 0)
			return -face_normal;
		else
			return face_normal;
	}
	//获得四面体的各种用于计算的梯度，其中vertex_index是遍历四面体的开始点
	void Nonsmoother::getGradient(int vertex_index, int hedron_index, HedronGradient &_hedron_gradient)
	{
		OneTetrahedron2 hedron;
		NonsmoothPoint point[4];
		NonsmoothPoint t, u, v;
		double volume;
		double edgelength[3][4];          /* the lengths of each of the edges of the tet */
		NonsmoothPoint edgegrad[3][4]; /* the gradient of each edge length wrt vtx1 */
		NonsmoothPoint facenormal[4];  /* the normals of each face of the tet */
		double facearea[4];               /* areas of the faces of the tet */
		NonsmoothPoint facegrad[4];    /* the gradient of each of the face areas wrt vtx1 */
		NonsmoothPoint volumegrad = NonsmoothPoint(0.0, 0.0, 0.0);
		NonsmoothPoint e1, e2;
		NonsmoothPoint term1, term2, diff;
		double top, bot;
		NonsmoothPoint top_gradient, bot_gradient, quotient_gradient;
		hedron = tetrahedrons2_[hedron_index];
		int i, j, k, l;
		double c;
		int edge_index = 0;
		//point[0] = points_[vertex_index];
		for (i = 0; i < 4; ++ i)
		{
			point[i] = points_[hedron.four_points_index[i]];
		}
		t = point[1] - point[0];
		v = point[2] - point[0];
		u = point[3] - point[0];
		//cout volume
		volume = ((t % v) | u) / 6.0;
		_hedron_gradient.volume = volume;
		for(i = 0; i < 4; ++ i)
		{
			//此处计算面法向，facenormal[i]是第i个点对面的面方向，此处法向方向的不对会对体积的梯度计算结果正负有影响
			//这个方法可以用自己写的getFacenormal来试一试
			j = (i + 1) & 3;
			if ((i & 1) == 0)
			{
				k = (i + 3) & 3;
				l = (i  + 2) & 3;
			}
			else
			{
				k = (i + 2) & 3;
				l = (i  + 3) & 3;
			}
			facenormal[i] = (point[k] - point[j]) % (point[l] - point[j]);
			if (i == 0)
			{
				_hedron_gradient.volumegrad = volumegrad = facenormal[i] /(- 6.0);
			}
			facearea[i] = sqrt(facenormal[i] | facenormal[i]) / 2.0;
			if (i == 0)
			{
				/* this face doesn't include vtx1, gradient is zero */
				facegrad[i] = NonsmoothPoint(0.0, 0.0, 0.0);
			}
			else
			{
				double factor = 1.0 / (4.0 * facearea[i]);
				switch(i)
				{
					/* compute the area of face 1 using u and v */
				case 1:
					e1 = u;
					e2 = v;
					break;
				case 2:
					e1 = t;
					e2 = u;
					break;
				case 3:
					e1 = v;
					e2 = t;
					break;
				}
				/* compute first term of gradient */
				c = e2 | (e1 - e2);
				term1 = NonsmoothPoint(c*e1[0], c*e1[1], c*e1[2]);

				/* compute the second term */
				c = e1 | (e1 - e2);
				term2 = NonsmoothPoint(c*e2[0], c*e2[1], c*e2[2]);

				/* now, combine the terms, scaled with the 1/4A */
				diff = term1 - term2;
				facegrad[i] = NonsmoothPoint(factor*diff[0], factor*diff[1], factor*diff[2]);
			}
			//MINSINE与VOLLENGTH需要计算边的长度和边的梯度
			if (quality_kinds_ == MINSINE2)
			{
				for (j = i + 1; j < 4; j++) 
				{
					/* e1 is edge from point i to point j */
					e1 = point[j] - point[i];
					edgelength[i][j] = sqrt(e1 | e1);

					/* also compute the gradient of the length of this edge */

					/* if vtx1 isn't one of the edge's endpoints, the gradent is zero */
					if (i != 0)
					{
						edgegrad[i][j] = NonsmoothPoint(0.0, 0.0, 0.0);
					}
					/* otherwise, it's in the negative direction of this edge,
					and scaled by edge length */
					else
					{ 
						double  factor = -1.0 / edgelength[i][j];
						edgegrad[i][j] = NonsmoothPoint(factor*e1[0], factor*e1[1], factor*e1[2]);
					}
				}
			}			
		}
		//#end compute the gradient of edgelength and area
		//#start compute quality gradient
		//misine quality gradient
		
		//MINSINE2需要计算sin以及sin的梯度
		if (quality_kinds_ == MINSINE2)
		{
			for(i = 0; i < 3; ++ i)
			{
				for (j = i + 1; j  < 4; ++ j)
				{
					k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
					l = 6 - i - j - k;
					/* compute the sine of this dihedral angle */
					//这边是不是有错误？？？根据公式应该是_hedron_gradient.sine[edge_index] = (3 * volume * edgelength[k][l]) / (2 * facearea[k] * facearea[l])
					//或者_hedron_gradient.sine[edge_index] = (3 * volume * edgelength[i][j]) / (2 * facearea[i] * facearea[j]);
					//相同地，top = volume * edgelength[i][j];是不是有错？？？
					//
					_hedron_gradient.sine[edge_index] = (3 * volume * edgelength[i][j]) / (2 * facearea[k] * facearea[l]);
					top = volume * edgelength[i][j];
					bot = facearea[k] * facearea[l];

					/* find gradient of top */
					top_gradient = gradproduct(volume, edgelength[i][j], volumegrad, edgegrad[i][j]);

					/* find gradient of bottom */
					bot_gradient = gradproduct(facearea[k], facearea[l], facegrad[k], facegrad[l]);

					/* now, find the gradient of the quotient */
					quotient_gradient = gradquotient(top, bot, top_gradient, bot_gradient);

					/* now scale with constant factor */
					_hedron_gradient.sinegrad[edge_index] = NonsmoothPoint(1.5*quotient_gradient[0], 1.5*quotient_gradient[1], 1.5*quotient_gradient[2]);
					edge_index++;
				}
			}
		}
	}
	void Nonsmoother::getActiveSet(int _vertex_index, double _worst_quality, NonsmoothPoint *_active_gradient, int & _active_num, HedronGradient * _gradients)
	{
		int incident_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		_active_num = 0;
		int i, j;

		for (i = 0; i < incident_num; ++ i)
		{
				for (j = 0; j < 6; j++)
				{
					/* is this close enough to the worst? */
					if (_gradients[i].sine[j] < (_worst_quality * ACTIVESETFACTOR))
					{
						/* get the actual gradient value */
						_active_gradient[_active_num] = _gradients[i].sinegrad[j];
						++ _active_num;
					}
				}
		}
	}
	void Nonsmoother::RecylefindBasis(NonsmoothPoint * _M, int _M_num, NonsmoothPoint * _S, int _S_num, NonsmoothPoint * _B, int &_B_num)
	{
		NonsmoothPoint p, q, r, s;                          /* vertices */
		NonsmoothPoint s1, t1, d1, d2;                      /* temporary vertices */
		NonsmoothPoint origin = NonsmoothPoint(0.0, 0.0, 0.0);
		NonsmoothPoint localS[MAXBASIS];             /* for passing to recursive calls */
		NonsmoothPoint localM[MAXBASIS];                 /* for passing to recursive colls */
		int localS_num, localM_num;
		while(true)
		{
			p = _M[_M_num-1];
			if (_M_num == 1)
			{
				/* and S has no elements */
				if (_S_num == 0)
				{
					/* then the single element in M must be the
					entire basis, just send back M */
					_B[_B_num] = _M[0];
					++ _B_num;
					break;
				}

				/* otherwise, because we assume the last element
				we just removed is not part of the basis, assign
				the basis to be the elements of S */
				for (int i = 0; i < _S_num; ++ i)
				{
					_B[_B_num] = _S[i];
					++ _B_num;
				}
			}
			else
			{
				/* make a new copy of M minus the last element */
				for(int i = 0; i < _M_num - 1; ++ i)
				{
					localM[i] = _M[i];
				}
				for(int i = 0; i < _S_num; ++ i)
				{
					localS[i] = _S[i];
				}
				localM_num = _M_num -1;
				localS_num = _S_num;
				continue;
			}
			/* if M has only one element */


			switch(_B_num)
			{
			case 1:
				q = _B[0];

				/* compute the vector from q to p */
				d1 = p - q;

				/* check the sign of the dot product. >=0 means p doesn't
				improve the basis */
				if ( (q | d1) >= 0)
				{
					/* this is a good B, send it back!*/
					continue;
				}
				break;
			case 2:
				/* fetch coordinates from the mesh */
				q = _B[0];
				r = _B[1];

				/* compute vector s from r to p */
				s1 = p - r;
				/* compute vector t from r to q */
				t1 = q - r;

				/* now a couple of cross products */
				d1 = s1 % t1;
				d2 = r % t1;

				/* can p improve the basis? */
				if ((d1 | d2) >= 0)
				{
					/* nope! send back B as is. */
					continue;
				}            
				break;
			case 3:
				/* fetch coordinates from the mesh */
				q = _B[0];
				r = _B[1];
				s = _B[2];

				/* does p improve the basis? */
				if (orient(p, q, r, s) * orient(origin, q, r, s) <= 0)
				{
					/* nope! send back B as is. */
					continue;
				}
				break;
			default:
				continue;
			}
			if ((_M_num == 1) || (_S_num == 3))
			{
				/* copy S into B */

				_B_num = 0;
				for (int i = 0; i < _S_num; ++ i)
				{
					_B[_B_num] = _S[i];
					++ _B_num;
				}
				_B[_B_num] = p;
				++ _B_num;
				continue;
			}
			else
			{
				/* create the new S */
				for (int i = 0; i < _S_num; ++ i)
				{
					localS[i] = _S[i];
				}
				localS_num = _S_num + 1;
				localS[_S_num] = p;

				/* create the new M, leaving off the last element */
				for (int i = 0; i < _M_num - 1; ++ i)
				{
					localM[i] = _M[i];
				}
				localM_num = _M_num - 1;

				/* find any basis points remaining in M */
				continue;
			}        
		}
	}
	void Nonsmoother::findBasis(NonsmoothPoint * _M, int _M_num, NonsmoothPoint * _S, int _S_num, NonsmoothPoint * _B, int &_B_num)
	{
		NonsmoothPoint p, q, r, s;                          /* vertices */
		NonsmoothPoint s1, t1, d1, d2;                      /* temporary vertices */
		NonsmoothPoint origin = NonsmoothPoint(0.0, 0.0, 0.0);
		NonsmoothPoint localS[MAXBASIS];             /* for passing to recursive calls */
		NonsmoothPoint localM[MAXBASIS];                 /* for passing to recursive colls */
		int localS_num, localM_num;
		p = _M[_M_num-1];
		int i;

		/* if M has only one element */
		if (_M_num == 1)
		{
			/* and S has no elements */
			if (_S_num == 0)
			{
				/* then the single element in M must be the
				entire basis, just send back M */
				_B[_B_num] = _M[0];
				++ _B_num;
				return;
			}

			/* otherwise, because we assume the last element
			we just removed is not part of the basis, assign
			the basis to be the elements of S */
			for (i = 0; i < _S_num; ++ i)
			{
				_B[_B_num] = _S[i];
				++ _B_num;
			}
		}
		else
		{
			/* make a new copy of M minus the last element */
			for(i = 0; i < _M_num - 1; ++ i)
			{
				localM[i] = _M[i];
			}
			for(i = 0; i < _S_num; ++ i)
			{
				localS[i] = _S[i];
			}
			localM_num = _M_num -1;
			localS_num = _S_num;
			findBasis(localM, localM_num, localS, localS_num, _B, _B_num);
		}
		switch(_B_num)
		{
		case 1:
			q = _B[0];

			/* compute the vector from q to p */
			d1 = p - q;

			/* check the sign of the dot product. >=0 means p doesn't
			improve the basis */
			if ( (q | d1) >= 0)
			{
				/* this is a good B, send it back!*/
				return;
			}
			break;
		case 2:
			/* fetch coordinates from the mesh */
			q = _B[0];
			r = _B[1];

			/* compute vector s from r to p */
			s1 = p - r;
			/* compute vector t from r to q */
			t1 = q - r;

			/* now a couple of cross products */
			d1 = s1 % t1;
			d2 = r % t1;

			/* can p improve the basis? */
			if ((d1 | d2) >= 0)
			{
				/* nope! send back B as is. */
				return;
			}            
			break;
		case 3:
			/* fetch coordinates from the mesh */
			q = _B[0];
			r = _B[1];
			s = _B[2];

			/* does p improve the basis? */
			if (orient(p, q, r, s) * orient(origin, q, r, s) <= 0)
			{
				/* nope! send back B as is. */
				return;
			}
			break;
		default:
			return;
		}
		if ((_M_num == 1) || (_S_num == 3))
		{
			/* copy S into B */

			_B_num = 0;
			for (int i = 0; i < _S_num; ++ i)
			{
				_B[_B_num] = _S[i];
				++ _B_num;
			}
			_B[_B_num] = p;
			++ _B_num;
			return;
		}
		else
		{
			/* create the new S */
			for (i = 0; i < _S_num; ++ i)
			{
				localS[i] = _S[i];
			}
			localS_num = _S_num + 1;
			localS[_S_num] = p;

			/* create the new M, leaving off the last element */
			for (i = 0; i < _M_num - 1; ++ i)
			{
				localM[i] = _M[i];
			}
			localM_num = _M_num - 1;

			/* find any basis points remaining in M */
			findBasis(localM, localM_num, localS, localS_num, _B, _B_num);
			return;
		}        
	}
	double Nonsmoother::getInitialAlpha(int _vertex_index, NonsmoothPoint _direciton, double _r, double _worst_quality, HedronGradient * _gradients)
	{
		double alpha = HUGEFLOAT;
		double newalpha;
		double rate;
		int i, j;
		int incident_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		if (_vertex_index == 44)
		{
			int x = 0;
			++ x;
		}
		for (i = 0; i < incident_num; ++ i)
		{
			for (j = 0; j < 6; j++)
			{
				/* if this function improves more slowly than any in the active set, 
				then it might end up as the objective. */
				rate = _direciton | _gradients[i].sinegrad[j];
				if (rate + RATEEPSILON < _r)
				{
					newalpha = (_gradients[i].sine[j] - _worst_quality) / (_r - rate);
					if (newalpha < alpha)
						alpha = newalpha;
				}						
			}
		}
		if (alpha < 0.0)
		{
			alpha = 0.0;
		}
		return alpha;
	}
	void Nonsmoother::findDirection(int _vertex_index, NonsmoothPoint * _active_gradient, int _active_num, NonsmoothPoint &_direciton)
	{
		NonsmoothPoint B[MAXBASIS];
		NonsmoothPoint TS[MAXBASIS];
		int S_num = 0;
		int B_num = 0;
		NonsmoothPoint p, q, r, pmq, s, t, s2, t2, sxt, rxt, sxr;
		double l, c, d;
		findBasis(_active_gradient, _active_num, TS, S_num, B, B_num);
		switch(B_num)
		{
		case 1:
			_direciton = B[0];
			return;
		case 2:
			p = B[0];
			q = B[1];
			pmq = p - q;
			l = sqrt(pmq | pmq);

			/* if these points are the same, just return one of them */
			if (l == 0.0)
			{
				_direciton = B[0];
				return;
			}

			c = (q | pmq) / (l * l);
			_direciton = c * pmq;
			_direciton = q - _direciton;
			return; 
		case 3:
			p = B[0];
			q = B[1];
			r = B[2];

			s = p - r;
			t = q - r;

			sxt = s % t;
			rxt = r % t;
			sxr = s % r;

			/* if any of these cross products is really tiny, give up
			and return the origin */
			if ((sqrt(sxt | sxt) < NEARESTMIN) || (sqrt(rxt | rxt) < NEARESTMIN) || (sqrt(sxr | sxr) < NEARESTMIN))
			{
				_direciton = NonsmoothPoint(0.0, 0.0, 0.0);
				return;
			}

			c = (sxt | rxt) / (sxt | sxt);
			d = (sxt | sxr) / (sxt | sxt);
			s2 = c *s;
			t2 = d * t;
			_direciton = r - s2 - t2;
			return;
		case 4:
			_direciton = NonsmoothPoint(0.0, 0.0, 0.0);
			return;
		}
	}
	void Nonsmoother::nonsmoothLineSearch(int _vertex_index, NonsmoothPoint _direction, double _r, double &_alpha, double _worst_quality)
	{
		int numiter = 0;              /* number of iterations */
		NonsmoothPoint origvertex;             /* to save the original vertex position */
		NonsmoothPoint offset;     /* the offset to move the vertex, alpha * d */
		double worstqual;             /* the current worst quality */
		double origworstqual;         /* the original worst quality */
		double thisqual;              /* the quality of the current tet */
		double oldworstqual;          /* the worst quality at the last step */
		NonsmoothPoint tetv[4];
		int hedron_index;
		OneTetrahedron2 hedron;
		origworstqual = oldworstqual = _worst_quality;
		origvertex = points_[_vertex_index];
		int hedron_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		int i, j;
		while ((_alpha > MINSTEPSIZE) && (numiter < MAXLINEITER))
		{
			offset = _alpha * _direction;
			points_[_vertex_index] = points_[_vertex_index] + offset;
			worstqual = HUGEFLOAT; 
			for (i = 0; i < hedron_num; i++)
			{
				hedron_index = incident_tetradrons2_[_vertex_index].incident_tetrahedrons[i];
				hedron = tetrahedrons2_[hedron_index];
				for (j = 0; j < 4; ++ j)
				{
					tetv[j] = points_[hedron.four_points_index[j]];
				}
				thisqual = nonsmoothquality_.getQuality(tetv, quality_kinds_);
				/* is this the worst quality we've seen? */
				if (thisqual < worstqual) worstqual = thisqual;
			}
			if ((oldworstqual > origworstqual) && (oldworstqual > worstqual))
			{
				/* use the previous step's alpha */
				_alpha = (_alpha) * 2;
				/* put vertex back where it started */
				points_[_vertex_index] = origvertex;
				return;
			}
			/* if we have succeeded in gaining 90% of the expected
			improvement, accept this initial step size */
			double factor = 0.9 * _alpha * _r;
			if ((worstqual - origworstqual) > factor)
			{
				/* put vertex back where it started */
				points_[_vertex_index] = origvertex;
				return;
			}
			/* put vertex back where it started */
			points_[_vertex_index] = origvertex;

			/* cut alpha down by half and try again */
			_alpha = (_alpha) / 2.0;

			/* save the worst quality from this step */
			oldworstqual = worstqual;
		}
		/* no positive alpha could be found that improved things... give up and return zero */
		if((worstqual - origworstqual) > 0)
			return;
		else
			_alpha = 0.0;
	}
	double Nonsmoother::getMinQuality(int _vertex_index)
	{
		double worst_qulity, this_quality;
		int hedron_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		NonsmoothPoint p[4];
		int hedron_index;
		OneTetrahedron2 hedron;
		int i, j;
		worst_qulity = HUGEFLOAT;
		for (i = 0; i < hedron_num; ++ i)
		{
			hedron_index = incident_tetradrons2_[_vertex_index].incident_tetrahedrons[i];
			hedron = tetrahedrons2_[hedron_index];
			for (j = 0; j < 4; ++ j)
			{
				p[j] = points_[hedron.four_points_index[j]];
			}
			this_quality = nonsmoothquality_.getQuality(p, quality_kinds_);
			if(this_quality < worst_qulity)
				worst_qulity = this_quality;
		}
		return worst_qulity;
	}
	bool Nonsmoother::laplacianSmooth(int _vertex_index, double & worst_quality, double thredholds)
	{
		int incident_points[MAXINCIDENTPOINT];
		int incident_hedron_num;
		int incident_point_num = 0;
		int point_index;
		int hedron_index;
		double old_quality, new_quality;
		int i, j, k;
		NonsmoothPoint original_point = points_[_vertex_index];
		NonsmoothPoint new_point = NonsmoothPoint(0.0, 0.0, 0.0);
		incident_hedron_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		for (i = 0; i < incident_hedron_num; ++ i)
		{
			hedron_index = incident_tetradrons2_[_vertex_index].incident_tetrahedrons[i];
			for (j = 0; j < 4; ++ j)
			{
				point_index = tetrahedrons2_[hedron_index].four_points_index[j];
				if (point_index != _vertex_index)
				{
					for (k = 0; k < incident_point_num; ++ k)
					{
						if (incident_points[k] == point_index)
							break;
					}
					if(k == incident_point_num)
					{
						incident_points[k] = point_index;
						++ incident_point_num;
					}
				}
			}
		}
		for (int i = 0; i < incident_point_num; ++ i)
		{
			new_point = new_point + points_[incident_points[i]];
		}
		new_point = new_point / incident_point_num;
		old_quality = getMinQuality(_vertex_index);
		points_[_vertex_index] = new_point;
		new_quality = getMinQuality(_vertex_index);
		if (new_quality > old_quality)
			return true;
		points_[_vertex_index] = original_point;
		return false;
	}
	bool Nonsmoother::NonsmoothVertice(int _vertex_index, double &_out_worst_quality)
	{
		int numiter = 0;       /* number of optimization iterations */
		double worstqual;      /* the numerical value of the worst quality function */
		double thisqual;       /* quality of the current tet */
		double oldworstqual;   /* the numerical value of the worst quality function at the last step */
		double oldest_qual;
		double improvement;    /* how much we've improved this step */
		double dlength;        /* length of d */
		double r;              /* expected rate of improvement */
		double newr;           /* temp r var */
		double alpha;          /* step size */
		double newalpha;       /* candidate new step size */
		double rate;           /* the rate a particular function will improve in the direction d */
		NonsmoothPoint change;      /* change we'll make to the point */
		NonsmoothPoint verts[4];        /* tetra points*/
		//SmoothVertex vinfo;    /* contains vertex type and a vector related to the vertex */

		//float allgrads[3];
		NonsmoothPoint d;           /* direction to move vertex */
		/* the gradients in the active set */
		NonsmoothPoint activegrads[MAXACTIVESET];
		NonsmoothPoint origpoint;       /* to save the original vertex location */
		NonsmoothPoint p = points_[_vertex_index];  /* the vertex to be altered */
		bool is_nonsmoothed = false;
		OneTetrahedron2 hedron;
		int hedron_index;
		int i;
		origpoint = p;
		worstqual = HUGEFLOAT; 
		int hedron_num = incident_tetradrons2_[_vertex_index].incident_hedron_num;
		HedronGradient * hedrongradient = new HedronGradient[hedron_num];

		for (i = 0; i < hedron_num; i++)
		{
			hedron_index = incident_tetradrons2_[_vertex_index].incident_tetrahedrons[i];
			getGradient(_vertex_index, hedron_index, hedrongradient[i]);
		}
		worstqual = getMinQuality(_vertex_index);
		
		if (worstqual < 0.0)
		{
			delete hedrongradient;
			return false;
		}
		oldest_qual = worstqual;
		int active_set_num = 0;
		getActiveSet(_vertex_index, worstqual, activegrads, active_set_num, hedrongradient);

		if(active_set_num == 0)
		{
			delete hedrongradient;
			return false;
		}
		//findDirection(_vertex_index, activegrads, active_set_num, d);
		NonsmoothPoint grad_point(0, 0, 0);
		for (i = 0; i < active_set_num; ++ i)
		{
			grad_point = grad_point + activegrads[i];
		}
		d = grad_point / active_set_num;
		dlength = sqrt(d | d);
		if(dlength < DEPSILON)
		{
			delete hedrongradient;
			return false;
		}
		//#begin nonsmooth
		do 
		{
			r = HUGEFLOAT;
			//这里计算r是做什么的？
			//这里d是不是应该先单位化呢？
			
			for (i = 0; i < active_set_num; ++ i)
			{
				newr = d | activegrads[i];
				if (newr <= 0.0)
				{
					delete hedrongradient;
					/* if we have already moved the vertex some, this is still a success */
					if (!(p == origpoint))
					{
						_out_worst_quality = worstqual;
						return true;
					}
					else
						return false;
				}
				if (newr < r)
				{
					r = newr;
				}
			}
			/* initialize alpha to the nearest intersection with another
			quality function */
			alpha = HUGEFLOAT;
			alpha = getInitialAlpha(_vertex_index, d, r, worstqual, hedrongradient);
			
			//如果alpha没有改变，我们要给alpha找一个设定初值？这里为何要使用体积？
			if (alpha == HUGEFLOAT)
			{
				for (i = 0; i < hedron_num; i++)
				{
					/* if moving in the direction d will decrease 
					this element's volume */					
					rate = d | hedrongradient[i].volumegrad;

					if (rate < 0.0)
					{
						newalpha = -hedrongradient[i].volume / (2.0 * rate);
						/* if this is smaller than the current step size,
						use it */
						if (newalpha < alpha)
							alpha = newalpha;
					}
				}
			}

			//寻找最佳的距离
			nonsmoothLineSearch(_vertex_index, d, r, alpha, worstqual);
			//用得到的结果替代当前点的位置然后判断是否有提升质量，有则进行下一次优化，没有则复回当前结果并
			change = alpha * d;

			if ((abs(change.point_coord[0]) > MAXCHANGE) || (abs(change.point_coord[1]) > MAXCHANGE) || (abs(change.point_coord[2]) > MAXCHANGE))
				break;
			
			//问题出在change中，在cow和dragon中change会达到10^16级别，这样会导致面积计算出问题导致出现质量为0的情况,求证的
			//结果是跟alpha关系很大，会达到10^16级别，d是正常的。所以，聚焦于getInitialAlpha，nonsmoothLineSearch两个函数，在出错
			//的时候，getinitialAlpha返回默认的最大值，所以，接下来看根据体积梯度来调整alpha是否正确
			//通过结果比对发现，用体积梯度调整的结果基本上使alpha变得更大甚至达到了10^35级别
			//如果getInitialAlpha会失效是正常的情况下，那么根据体积梯度调整alpha的方法应该有问题吧
			//重点查看getInitialAlpha和梯度调整法的算法
			//体积的计算大多正常，可是有的竟然达到10^34级别，体积按道理应该<100比较正常
			//经过分析，是因为传入的点具有10^16级别的数值导致体积计算错误，
			//所以体积的计算方法是没有错的，只是传入的数据本身有问题
			//那么问题究竟出在哪里呢？
			//有待解决，临时的解决方法是限定一个MAXCHANGE(这里为100)，使过大的改变不进行处理
			//另一个想法是change= d，发现结果和不处理基本一样。

			p = p + change;			
			points_[_vertex_index] = p;
			oldworstqual = worstqual;
			worstqual = getMinQuality(_vertex_index);
			
			//for test
			improvement = worstqual - oldworstqual;

			//if (improvement > MINSMOOTHITERIMPROVE)
			//	is_nonsmoothed = true;
			//assert(improvement >= 0);
			if (improvement <= 0)
			{
				p = p - change;
				points_[_vertex_index] = p;
				worstqual = oldworstqual;
				break;
			}

			for(i = 0; i < hedron_num; ++ i)
			{
				hedron_index = incident_tetradrons2_[_vertex_index].incident_tetrahedrons[i];
				getGradient(_vertex_index, hedron_index, hedrongradient[i]);
			}

			active_set_num = 0;
			getActiveSet(_vertex_index, worstqual, activegrads, active_set_num, hedrongradient);
			if(active_set_num == 0)
			{
				delete hedrongradient;
				return false;
			}
			//findDirection(_vertex_index, activegrads, active_set_num, d);
			
			for (i = 0; i < active_set_num; ++ i)
			{
				grad_point = grad_point + activegrads[i];
			}
			d = grad_point / active_set_num;
			dlength = sqrt(d | d);
			++ numiter;
		} while ((dlength > DEPSILON) && 
			(numiter < MAXSMOOTHITER) &&
			(improvement > MINSMOOTHITERIMPROVE));

		_out_worst_quality = worstqual;
		delete hedrongradient;

		if (worstqual - oldest_qual > 0)
			return true;
		else
			return false;
	}
	bool Nonsmoother::NonsmoothVertice(double p[3], double incidenttet[][9], int tetctn, double newp[3])
	{
		std::set<NonsmoothPoint> points;
		std::set<NonsmoothPoint>::iterator sp_iter;
		for (int i = 0; i < tetctn; ++ i)
		{
			for (int j = 0; j < 3; ++ j)
			{
				NonsmoothPoint temp_point(incidenttet[i][3 * j], incidenttet[i][3 * j + 1], incidenttet[i][3 * j + 2]);
				points.insert(temp_point);
			}
		}
		points_ = new NonsmoothPoint[1 + points.size()];
		incident_tetradrons2_ = new IncidentTetrahedrons2[1];
		tetrahedrons2_ = new OneTetrahedron2[points.size()];
		//初始化所有的顶点
		points_[0] = NonsmoothPoint(p);
		int point_index = 1;
		for (sp_iter = points.begin(); sp_iter != points.end(); ++ sp_iter)
		{
			points_[point_index] = * sp_iter;
			++ point_index;
		}
		vertex_num_ = point_index;
		hedron_num_ = tetctn;
		incident_tetradrons2_[0].incident_hedron_num = tetctn;
		//初始化所有四面体
		for (int i = 0; i < tetctn; ++ i)
		{
			OneTetrahedron2 hedron;
			hedron.four_points_index[0] = 0;
			int four_index = 1;
			for (int j = 0; j < 3; ++ j)
			{
				NonsmoothPoint temp_point(incidenttet[i][3 * j], incidenttet[i][3 * j + 1], incidenttet[i][3 * j + 2]);
				for (int k = 1; k < vertex_num_; ++ k)
					if (points_[k] == temp_point)
					{
						hedron.four_points_index[four_index] = k;
						++ four_index;
						break;
					}
				}
			tetrahedrons2_[i] = hedron;
			//初始化该点的邻接四面体的index
			incident_tetradrons2_[0].incident_tetrahedrons[i] = i;
			}
		double temp_quality;
		bool suc = NonsmoothVertice(0, temp_quality);
		for ( int i = 0; i < 3; ++ i)
		{
			newp[i] = points_[0].point_coord[i];
		}
		return suc;
	}
		
	bool Nonsmoother::NonsmoothVertex(double thredhold)
	{
		double quality_after = 0;
		double quality_before = 0;
		double worst_quality = 0;
		smooth_time_ = 0;
		do 
		{
			quality_before = quality_after;
			++ smooth_time_;
			worst_quality = 1.0e-5;
			
			for (int i = 0; i < vertex_num_; ++ i)
			{
				if (boundaryChange_ || (!points_[i].boundary))
				{
					switch (smoothing_type_)
					{
					case 0x01:laplacianSmooth(i, worst_quality, MINSMOOTHITERIMPROVE);break;
					case 0x10:NonsmoothVertice(i, worst_quality);break;
					case 0x11:{
						if(!laplacianSmooth(i, worst_quality, MINSMOOTHITERIMPROVE))
							NonsmoothVertice(i, worst_quality);
							  }
					}
				}
		}			
		} while ( (smooth_time_ < SMOOTHTIME) );
		//这里循环没有限制(abs(change.point_coord[0]) > MAXCHANGE)，导致rand2会出错，质量为0，原因不详，应该是change
		//的值限制不够。
		if (new_quality_ - old_quality_ > thredhold)
			return true;
		else
			return false;
	}
	void Nonsmoother::saveData(std::string _file_name)
	{
		std::ofstream outfile;
		int i, j;
		std::string file_name = "F:\\SaveModel\\" + _file_name + ".txt";
		outfile.open(file_name.data());
		outfile<<vertex_num_<<std::endl;
		outfile<<hedron_num_<<std::endl;
		/*for (i = 0; i < vertex_num_; ++ i)
		{
			outfile<<points_[i].point_coord[0]<<" "<<points_[i].point_coord[1]<<" "<<points_[i].point_coord[2];
			if (points_[i].boundary)
				outfile<<" 1"<<std::endl;
			else
				outfile<<" 0"<<std::endl;
		}*/

		TetraMesh::PointIter p_iter;
		TetraMesh::PointHandle ph;
		int index;
		TetraMesh::Point point;
		for (p_iter = tmesh_->points_begin(); p_iter != tmesh_->points_end(); ++ p_iter)
		{
			ph = p_iter.handle();
			index = ph.idx();
			point = tmesh_->point(ph);
			outfile<<point[0]<<" "<<point[1]<<" "<<point[2];
			if (tmesh_->is_boundary(ph))
				outfile<<" 1"<<std::endl;
			else
				outfile<<" 0"<<std::endl;
		}


		for (i = 0; i < hedron_num_; ++ i)
		{
			for (j = 0; j < 4; ++ j)
			{
				outfile<<tetrahedrons2_[i].four_points_index[j]<<" ";
			}
			if(tetrahedrons2_[i].boundary)
				outfile<<"1"<<std::endl;
			else
				outfile<<"0"<<std::endl;
		}
		outfile.close();
	}
	void Nonsmoother::readData(std::string _file_name)
	{
		std::ifstream infile;
		int i, j, k, l;
		double point_coord;
		std::string file_name = "F:\\saveModel\\" + _file_name + ".txt";
		infile.open(file_name.data());
		infile>>vertex_num_;
		infile>>hedron_num_;
		points_ = new NonsmoothPoint[vertex_num_];
		old_points_ = new NonsmoothPoint[vertex_num_];
		tetrahedrons2_ = new OneTetrahedron2[hedron_num_];
		incident_tetradrons2_ = new IncidentTetrahedrons2[vertex_num_];
		for (i = 0; i < vertex_num_; ++ i)
		{
			for (j = 0; j < 3; ++ j)	
			{
				infile>>point_coord;
				points_[i].point_coord[j] = old_points_[i].point_coord[j] = point_coord;
			}
			infile>>k;
			if (k == 1)	
			{
				points_[i].boundary = true;
				old_points_[i].boundary = true;
			}
			else
			{
				points_[i].boundary = false;
				old_points_[i].boundary = false;
			}
			incident_tetradrons2_[i].incident_hedron_num = 0;
		}
		for (i = 0; i < hedron_num_; ++ i)
		{
			for (j = 0; j < 4; ++ j)
			{
				infile>>k;
				tetrahedrons2_[i].four_points_index[j] = k;
				l = incident_tetradrons2_[k].incident_hedron_num;
				incident_tetradrons2_[k].incident_tetrahedrons[l] = i;
				++ incident_tetradrons2_[k].incident_hedron_num;
			}
			infile>>k;
			if (k == 1)
				tetrahedrons2_[i].boundary = true;
			else
				tetrahedrons2_[i].boundary = false;
		}
		infile.close();
	}
	void Nonsmoother::setValue(int _vertex_num, int _hedron_num, double * _points, int * _hedron, int * _incident_hedron, int * _incident_hedron_num, int * _incident_point, int * _incident_point_num)
	{
		vertex_num_ = _vertex_num;
		hedron_num_ = _hedron_num;
		points_ = new NonsmoothPoint[vertex_num_];
		old_points_ = new NonsmoothPoint[vertex_num_];
		tetrahedrons2_ = new OneTetrahedron2[hedron_num_];
		incident_tetradrons2_ = new IncidentTetrahedrons2[vertex_num_];
		int i, j;
		for (i = 0; i < hedron_num_; ++ i)
		{
			for(j = 0; j < 4; ++ j)
			{
				tetrahedrons2_[i].four_points_index[j] = _hedron[4 * i + j];
			}
		}
		for (i = 0; i < vertex_num_; ++ i)
		{
			for (j = 0; j < 3; ++ j)
			{
				points_[i].point_coord[j] = _points[3 * i + j];
			}
			incident_tetradrons2_[i].incident_hedron_num = _incident_hedron_num[i];
			for (j = 0; j < _incident_hedron_num[i]; ++ j)
			{
				incident_tetradrons2_[i].incident_tetrahedrons[j] = _incident_hedron[300 * i + j];
			}
		}
	}
	void Nonsmoother::pointGroup()
	{
		int color_vector[50];
		point_color_ = new int[vertex_num_];
		group_num_ = 0;
		int vertex_num;
		std::set<int> incident_point_num;
		std::set<int>::iterator si_iter;
		for(int i = 0; i < vertex_num_; ++ i)
		{
			//分组颜色初始化
			point_color_[i] = -1;
		}
		//对点进行分组，为并行使用
		for(int i = 0; i < vertex_num_; ++ i)
		{
			memset(color_vector, 0, sizeof(int)*group_num_);
			//根据邻接四面体得到邻接点集合
			int incident_num = incident_tetradrons2_[i].incident_hedron_num;
			incident_point_num.clear();
			for (int j = 0; j < incident_num; ++ j)
			{
				int hedron_num = incident_tetradrons2_[i].incident_tetrahedrons[j];
				for (int k = 0; k < 4; ++ k)
				{
					incident_point_num.insert(tetrahedrons2_[hedron_num].four_points_index[k]);
				}
			}
			for (si_iter = incident_point_num.begin(); si_iter != incident_point_num.end(); ++ si_iter)
			{
				vertex_num = *si_iter;
				if (point_color_[vertex_num] != -1)
				{
					color_vector[point_color_[vertex_num]] = 1;
				}
			}
			for (int j = 0; j < group_num_; ++j)
			{
				if (!color_vector[j])
				{
					point_color_[i] = j;
					break;
				}
			}
			// if exist color all used then add a new color
			if (point_color_[i] == -1)
			{
				point_color_[i] = group_num_ ++;
			}
		}
	}
}