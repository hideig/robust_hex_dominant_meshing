#include <TopologyProcess/EdgeContract.h>
#include <TopologyProcess/TopologyFlip.h>

namespace EdgeContractProcess
{

	void EdgeContract::initialize()
	{
		if (!mesh)
			return;

		// build edge
		char buf[128];
		VolumeMesh::PointHandle p[2];
		VolumeMesh::TetraMesh::HedronIter h_it;
		VolumeMesh::TetraMesh::HedronHalfEdgeIter hhe_it;
		edgeVec.clear();
		_ephc.clear();
		for (h_it = mesh->hedrons_begin(); h_it != mesh->hedrons_end(); ++ h_it)
		{
			for (hhe_it = mesh->hedron_half_edge_iter(h_it); hhe_it; ++hhe_it)
			{
				p[0] = mesh->point_handle(mesh->from_vertex_handle(hhe_it.handle()));
				p[1] = mesh->point_handle(mesh->to_vertex_handle(hhe_it.handle()));
				std::sort(p, p + 2);
				sprintf_s(buf,"%d %d",p[0],p[1]);
				if (_ephc.find(buf) == _ephc.end())
				{
					_ephc.insert(buf);
					EdgeTH edge(hhe_it.handle());
					edgeVec.push_back(edge);
				}
			}
		}
	}

	void EdgeContract::edgeContractData()
	{
		initialize();

		// open data saving file
		bool isBoundary;
		std::string fileName, fileName1, fileName2, fileName3;
		fileName = "F:\\ec_data\\ec_" + entityName;
		fileName1 = "F:\\ec_data\\ec_success";
		fileName2 = "F:\\ec_data\\ec_fail";
		fileName3 = "F:\\ec_data\\ec_time_consumption";
		outfile.open(fileName.c_str());
		outfile1.open(fileName1.c_str(), std::ofstream::app);
		outfile2.open(fileName2.c_str(), std::ofstream::app);
		outfile_time.open(fileName3.c_str(), std::ofstream::app);
		
		outfile << "                                                         \n";
		// variables for data saving
		int ec_flag;
		int es_tn, pfs_tn, pts_tn, ps_ec_tn, er_tn;
		int ec_num = 0, ec_fail_num = 0;
		int location;    // record the worst tetra is round edge(0) or endpoints(1)
		double qVal_best, qVal_worst, qVal_ec_best, qVal_ec_worst, qvb_, qvw_, q_improve, qVal_min_sine, qVal_square_root;
		// side ratio incident to endpoints, contract edge and all of them
		double sr_ep_small, sr_ep_large, sr_ep_aver, sr_ec_small, sr_ec_large, sr_ec_aver, sr_small, sr_large, sr_aver;
		VolumeMesh::TetraMesh::HalfEdgeHandle heh_;
		double min_dihedral_sine;
		std::vector<double> aroundJacVal;

		VolumeMesh::PointHandle ph_new;
		std::vector<VolumeMesh::TetraMesh::HedronHandle> pointStar;
		std::vector<VolumeMesh::TetraMesh::HedronHandle> edgeStar;

		c_ec = 0;
		c_ec_undo = 0;
		c_ec_extra = 0;
		c_ec_quality = 0;

		if (!QueryPerformanceFrequency(&li_start))
			cout << "QueryPerformanceFrequency failed\n";
		PCFreq = double(li_start.QuadPart);
		QueryPerformanceCounter(&li_start);
		for (unsigned int i = 0; i < edgeVec.size(); i++)
		{
			edgeVec[i].set_suc_flag(false);
			isBoundary = false;
			edgeStar.clear();
			mesh->edge_star(edgeVec[i].half_edge_handle(), edgeStar);
			// check if the edge is an boundary edge
			for (unsigned int j = 0; j < edgeStar.size(); j++)
			{
				if (mesh->is_boundary(edgeStar[j]))
				{
					isBoundary = true;
					break;
				}
			}
			if (isBoundary)
				continue;


			//start_t = clock();
			QueryPerformanceCounter(&li_start_t);
			// data saving before edge contractint
			++ ec_num;
			es_tn = edgeStar.size();
			best_worst_quality(edgeStar, qVal_best, qVal_worst, Q_VOL_LEN);
			best_worst_quality(edgeStar, qvb_, qvw_, Q_MIN_SIN);
			qVal_min_sine = 1.0;
			if (qVal_min_sine > qvw_)
				qVal_min_sine = qvw_;
			best_worst_quality(edgeStar, qvb_, qvw_, Q_SQUARE_ROOT);
			qVal_square_root = 1.0;
			if (qVal_square_root > qvw_)
				qVal_square_root = qvw_;

			// quality calculate time consumption
			QueryPerformanceCounter(&li_start_q);
			best_worst_quality(edgeStar, qvb_, qvw_, Q_MIN_SIN);
			QueryPerformanceCounter(&li_finish_q);
			c_ec_quality += double(li_start_q.QuadPart - li_finish_q.QuadPart);

			// side ratio of edge star
			sr_ec_small = sr_small = 10e6;
			sr_ec_large = sr_large = sr_ec_aver = sr_aver = 0;
			for (unsigned int k = 0; k < edgeStar.size(); k ++)
			{
				double t = side_ratio(edgeStar[k]);
				sr_aver += t;
				sr_ec_aver += t;
				if (sr_ec_large < t)
					sr_ec_large = t;
				if (sr_ec_small > t)
					sr_ec_small = t;
				if (sr_large < t)
					sr_large = t;
				if (sr_small > t)
					sr_small = t;
			}
			sr_ec_aver /= edgeStar.size();

			// get tetras containing edge's from_vertex
			location = 0;
			pointStar.clear();
			mesh->vertex_star(mesh->from_vertex_handle(edgeVec[i].half_edge_handle()), pointStar);
			pfs_tn = pointStar.size() - edgeStar.size();
			best_worst_quality(pointStar, qvb_, qvw_, Q_VOL_LEN);
			if (qVal_best < qvb_)
				qVal_best = qvb_;
			if (qVal_worst > qvw_)
			{
				qVal_worst = qvw_;
				location = 1;
			}
			best_worst_quality(pointStar, qvb_, qvw_, Q_MIN_SIN);
			if (qVal_min_sine > qvw_)
				qVal_min_sine = qvw_;
			best_worst_quality(pointStar, qvb_, qvw_, Q_SQUARE_ROOT);
			if (qVal_square_root > qvw_)
				qVal_square_root = qvw_;

			// quality calculate time consumption
			QueryPerformanceCounter(&li_start_q);
			best_worst_quality(pointStar, qvb_, qvw_, Q_MIN_SIN);
			QueryPerformanceCounter(&li_finish_q);
			c_ec_quality += double(li_start_q.QuadPart - li_finish_q.QuadPart);

			// side ratio of endpoint star
			sr_ep_aver = sr_ep_large = 0;
			sr_ep_small = 10e6;
			for (unsigned int k = 0; k < pointStar.size(); k++)
			{
				if (std::find(edgeStar.begin(), edgeStar.end(), pointStar[k]) != edgeStar.end())
					continue;

				double t = side_ratio(pointStar[k]);
				sr_ep_aver += t;
				sr_aver += t;
				if (sr_ep_large < t)
					sr_ep_large = t;
				if (sr_ep_small > t)
					sr_ep_small = t;
				if (sr_large < t)
					sr_large = t;
				if (sr_small > t)
					sr_small = t;
			}

			// get tetras containing edge's to_vertex
			pointStar.clear();
			mesh->vertex_star(mesh->to_vertex_handle(edgeVec[i].half_edge_handle()), pointStar);
			pts_tn = pointStar.size() - edgeStar.size();
			best_worst_quality(pointStar, qvb_, qvw_, Q_VOL_LEN);
			if (qVal_best < qvb_)
				qVal_best = qvb_;
			if (qVal_worst > qvw_)
			{
				qVal_worst = qvw_;
				location = 1;
			}
			best_worst_quality(pointStar, qvb_, qvw_, Q_MIN_SIN);
			if (qVal_min_sine > qvw_)
				qVal_min_sine = qvw_;
			best_worst_quality(pointStar, qvb_, qvw_, Q_SQUARE_ROOT);
			if (qVal_square_root > qvw_)
				qVal_square_root = qvw_;

			er_tn = es_tn + pfs_tn + pts_tn;

			// quality calculate time consumption
			QueryPerformanceCounter(&li_start_q);
			best_worst_quality(pointStar, qvb_, qvw_, Q_MIN_SIN);
			QueryPerformanceCounter(&li_finish_q);
			c_ec_quality += double(li_start_q.QuadPart - li_finish_q.QuadPart);

			// side ratio of endpoint star
			for (unsigned int k = 0; k < pointStar.size(); k++)
			{
				if (std::find(edgeStar.begin(), edgeStar.end(), pointStar[k]) != edgeStar.end())
					continue;

				double t = side_ratio(pointStar[k]);
				sr_ep_aver += t;
				sr_aver += t;
				if (sr_ep_large < t)
					sr_ep_large = t;
				if (sr_ep_small > t)
					sr_ep_small = t;
				if (sr_large < t)
					sr_large = t;
				if (sr_small > t)
					sr_small = t;
			}
			sr_aver /= er_tn;
			sr_ep_aver /= pfs_tn + pts_tn;

			// minimum dihedral_sine of edgestar
			min_dihedral_sine = 1;
			VolumeMesh::TetraMesh::HedronHalfEdgeIter hhe_it;
			double dihedral_sine_t;
			for (unsigned int k = 0; k < edgeStar.size(); k++)
			{
				VolumeMesh::VertexHandle pf,pt, pf_curr, pt_curr;
				pf = mesh->from_vertex_handle(edgeVec[i].half_edge_handle());
				pt = mesh->to_vertex_handle(edgeVec[i].half_edge_handle());
				for (hhe_it = mesh->hedron_half_edge_iter(edgeStar[k]); hhe_it; ++hhe_it)
				{
					pf_curr = mesh->from_vertex_handle(hhe_it.handle());
					pt_curr = mesh->to_vertex_handle(hhe_it.handle());
					if ((pf == pf_curr && pt == pt_curr) || (pf == pt_curr && pt == pf_curr))
					{
						dihedral_sine_t = dihedral_sine(hhe_it.handle());
						break;
					}
				}
				if (min_dihedral_sine > dihedral_sine_t)
				{
					min_dihedral_sine = dihedral_sine_t;
				}
			}

			//收缩边一周的四面体网格的雅克比的值
			aroundJacVal.clear();
			jac_val_around_edge(edgeVec[i].half_edge_handle(),aroundJacVal);

			QueryPerformanceCounter(&li_finish_t);
			c_ec_extra += double(li_finish_t.QuadPart - li_start_t.QuadPart);
 
			//---------------------------------------------------------------------------------------------------
			// do edge contract
			mesh->edge_contract(edgeVec[i].half_edge_handle(), ph_new);
			//----------------------------------------------------------------------------------------------------

			// addition
			QueryPerformanceCounter(&li_start_t);
			pointStar.clear();
			mesh->point_star(ph_new, pointStar);

			// data saving after operation
			ps_ec_tn = pointStar.size();
			best_worst_quality(pointStar, qVal_ec_best, qVal_ec_worst, Q_VOL_LEN);

			// quality calculate time consumption
			QueryPerformanceCounter(&li_start_q);
			best_worst_quality(pointStar, qvb_, qvw_, Q_MIN_SIN);
			QueryPerformanceCounter(&li_finish_q);
			c_ec_quality += double(li_start_q.QuadPart - li_finish_q.QuadPart);

			q_improve = qVal_ec_worst - qVal_worst;

			// set ec flag
			if (qVal_worst > qVal_ec_worst)
			{
				++ ec_fail_num;
				ec_flag = 0;
				edgeVec[i].set_suc_flag(false);
			}
			else
			{
				ec_flag = 1;
				edgeVec[i].set_suc_flag(true);
			}

			// output data
			outfile << ec_flag << "\t" << pfs_tn << "\t" << pts_tn << "\t" << es_tn << "\t" << er_tn << "\t"
				    << qVal_best << "\t" << qVal_worst << "\t" << ps_ec_tn << "\t"
					<< qVal_ec_best << "\t" << qVal_ec_worst << "\t" << location << "\t" << q_improve << "\t"
					<< sr_large    << "\t" << sr_small    << "\t" << sr_aver    << "\t"
					<< sr_ec_large << "\t" << sr_ec_small << "\t" << sr_ec_aver << "\t"
					<< sr_ep_large << "\t" << sr_ep_small << "\t" << sr_ep_aver << "\t"
					<< qVal_min_sine << "\t"  << qVal_square_root << "\t";

			for (unsigned int i = 0; i < aroundJacVal.size(); i++)
			{
				outfile << aroundJacVal[i] << "\t";
			}
			outfile << "\n";

			if (ec_flag)
			{
				outfile1 << ec_flag << "\t" << pfs_tn << "\t" << pts_tn << "\t" << es_tn << "\t" << er_tn << "\t"
				         << qVal_best << "\t" << qVal_worst << "\t" << ps_ec_tn << "\t"
						 << qVal_ec_best << "\t" << qVal_ec_worst << "\t" << location << "\t" << q_improve << "\t"
						 << sr_large    << "\t" << sr_small    << "\t" << sr_aver    << "\t"
						 << sr_ec_large << "\t" << sr_ec_small << "\t" << sr_ec_aver << "\t"
						 << sr_ep_large << "\t" << sr_ep_small << "\t" << sr_ep_aver << "\t"
						 << qVal_min_sine << "\t"  << qVal_square_root << "\t";

				for (unsigned int i = 0; i < aroundJacVal.size(); i++)
				{
					outfile1 << aroundJacVal[i] << "\t";
				}
				outfile1 << "\n";
			}
			else
			{
				outfile2 << ec_flag << "\t" << pfs_tn << "\t" << pts_tn << "\t" << es_tn << "\t" << er_tn << "\t"
				         << qVal_best << "\t" << qVal_worst << "\t" << ps_ec_tn << "\t"
				         << qVal_ec_best << "\t" << qVal_ec_worst << "\t" << location << "\t" << q_improve << "\t"
				         << sr_large    << "\t" << sr_small    << "\t" << sr_aver    << "\t"
				         << sr_ec_large << "\t" << sr_ec_small << "\t" << sr_ec_aver << "\t"
						 << sr_ep_large << "\t" << sr_ep_small << "\t" << sr_ep_aver << "\t"
						 << qVal_min_sine << "\t"  << qVal_square_root << "\t";

				for (unsigned int i = 0; i < aroundJacVal.size(); i++)
				{
					outfile2 << aroundJacVal[i] << "\t";
				}
				outfile2 << "\n";
			}

			QueryPerformanceCounter(&li_finish_t);
			c_ec_extra += double(li_finish_t.QuadPart - li_start_t.QuadPart);

			QueryPerformanceCounter(&li_start_t);
			mesh->recover_edge_contract();
			QueryPerformanceCounter(&li_finish_t);
			c_ec_extra += double(li_finish_t.QuadPart - li_start_t.QuadPart);
		}

		QueryPerformanceCounter(&li_finish);
		c_ec = double(li_finish.QuadPart - li_start.QuadPart);

		outfile.seekp(ios::beg);
		outfile << ec_num << " " << ec_num - ec_fail_num ;
		outfile_time << entityName << "\t" 
			         << (c_ec+c_ec_quality-c_ec_extra)/PCFreq << "\t" 
		             << (c_ec+c_ec_quality-c_ec_extra)/PCFreq/edgeVec.size() << "\n";

		outfile.close();
		outfile1.close();
		outfile2.close();
		outfile_time.close();
	}


	//---------------------------------hedron quality calculate------------------------------------//

	void EdgeContract::best_worst_quality(std::vector<VolumeMesh::TetraMesh::HedronHandle> tetraVec, 
		                                  double &qVal_best, double &qVal_worst, int Q_Type)
	{
		if (!tetraVec.size())
			return;

		qVal_best = 0;
		qVal_worst = 1;
		double qVal;
		VolumeMesh::Point p[4];
		VolumeMesh::TetraMesh::HedronVertexIter hv_it;
		for (unsigned int i = 0; i < tetraVec.size(); i++)
		{
			hv_it = mesh->hedron_vertex_iter(tetraVec[i]);
			p[0] = mesh->point(hv_it.handle());
			++ hv_it;
			p[1] = mesh->point(hv_it.handle());
			++ hv_it;
			p[2] = mesh->point(hv_it.handle());
			++ hv_it;
			p[3] = mesh->point(hv_it.handle());

			qVal = hedron_quality(p, Q_Type);
			if (qVal > qVal_best)
				qVal_best = qVal;
			if (qVal < qVal_worst)
				qVal_worst = qVal;
		}
	}

	double EdgeContract::hedron_quality(VolumeMesh::Point p[4], int Q_Type/* = Q_VOL_LEN*/)
	{
		double qVal;
		if (Q_Type == Q_VOL_LEN)
		{
			qVal = volume_length(p);
		}
		else if (Q_Type == Q_MIN_SIN)
		{
			//qVal = minimum_sine(p);
			qVal = minimum_sine_new(p);
		}
		else if (Q_Type == Q_JAC_VAL)
		{
			qVal = minimun_jacobian(p);
		}
		else if (Q_Type == Q_SQUARE_ROOT)
		{
			qVal = square_root(p);
		}
		return qVal;
	}


	//Volume-length measure, the best quality value is 1
	double EdgeContract::volume_length(VolumeMesh::Point v[4])
	{
		VolumeMesh::TetraMesh::Normal e[6];
		double volume;
		double l_rms = 0;
		e[0] = v[1] - v[0];
		e[1] = v[2] - v[0];
		e[2] = v[3] - v[0];
		e[3] = v[1] - v[2];
		e[4] = v[1] - v[3];
		e[5] = v[3] - v[2];

		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				l_rms += e[i][j] * e[i][j];
			}
		}
		l_rms /= 6.0;
		l_rms = pow(l_rms, 0.5);

		VolumeMesh::TetraMesh::Normal temp;
		temp = e[0] % e[1];
		volume = ((e[0] % e[1]) | e[2]) / 6;

		return (6 * pow(2, 0.5) * volume / pow(l_rms, 3));
	}

	// minimum sine measure
	double EdgeContract::minimum_sine(VolumeMesh::Point v[4])
	{
		double min_sine = 1, temp_sine;
		for (int i = 0; i < 3; i ++)
		{
			for (int j = i + 1; j < 4; j ++)
			{
				temp_sine = dihedral_sine(v, i, j);
				if (temp_sine < min_sine)
					min_sine = temp_sine;
			}
		}
		return min_sine;
	}

	double EdgeContract::minimum_sine_new(VolumeMesh::Point p[4])
	{
		double area[4];
		VolumeMesh::TetraMesh::Normal t, u, v;
		t = p[1] - p[0];
		u = p[2] - p[0];
		v = p[3] - p[0];
		area[0] = ((u - v) % (t - v)).norm()/2;
		area[1] = (u % v).norm()/2;
		area[2] = (v % t).norm()/2;
		area[3] = (t % u).norm()/2;
		double volume;
		volume = fabs((u%v)|t)/6;
		double minSine = 10e6;
		double temp;
		for (int i = 0; i < 3; i++)
		{
			for (int j = i+1; j < 4; j ++)
			{
				temp = (p[j] - p[i]).norm()/area[i]/area[j];
				if (minSine > temp)
					minSine = temp;
			}
		}
		minSine = minSine * 1.5 * volume;
		return minSine;
	}

	double EdgeContract::square_root(VolumeMesh::Point p[4])
	{
		double area[4];
		//VolumeMesh::TetraMesh::Normal t, u, v;
		VolumeMesh::Vec3d t, u, v;
		t = p[1] - p[0];
		u = p[2] - p[0];
		v = p[3] - p[0];
		area[0] = ((u - v) % (t - v)).norm()/2;
		area[1] = (u % v).norm()/2;
		area[2] = (v % t).norm()/2;
		area[3] = (t % u).norm()/2;
		double volume;
		volume = fabs((u%v)|t)/6;

		VolumeMesh::TetraMesh::Normal z[3];
		double zp[3];
		z[0] = u % v;
		z[1] = v % t;
		z[2] = t % u;
		zp[0] = pow(t.norm(),2);
		zp[1] = pow(u.norm(),2);
		zp[2] = pow(v.norm(),2);
		for (int i = 0; i < 3; i ++)
		{
			z[i][0] *= zp[i];
			z[i][1] *= zp[i];
			z[i][2] *= zp[i];
		}

		double Z;
		Z = (z[0] + z[1] + z[2]).norm();

		double sr;
		sr = 6 * pow(3,0.5) * volume / pow(Z*(area[0]+area[1]+area[2]+area[3]),0.5);
		return sr;
	}

	// calculate dihedral sine
	/** param[v]: points of tetrahedron
	param[k, l]: indexes of vertex of the shared edge
	**/
	double EdgeContract::dihedral_sine(VolumeMesh::Point v[4], int k, int l)
	{
		VolumeMesh::TetraMesh::Normal n[4];
		n[0] = (v[1] - v[3]) % (v[2] - v[3]);
		n[1] = (v[3] - v[0]) % (v[2] - v[0]);
		n[2] = (v[1] - v[0]) % (v[3] - v[0]);
		n[3] = (v[2] - v[0]) % (v[1] - v[0]);

		// get face index
		int fk, fl;
		if (fabs((k - l) * 1.0) == 2.0)
		{
			fk = (k + 1) % 4;
			fl = (l + 1) % 4;
		}
		else
		{
			fk = (k + 2) % 4;
			fl = (l + 2) % 4;
		}

		double sine, bestSine;
		bestSine = pow(3.0, 0.5) / 2.0;  // sin(60°)
		sine = pow(1 - pow((n[fk] | n[fl]) / n[fk].norm() / n[fl].norm(), 2.0), 0.5);

		if ((sine / bestSine) > 1)
		{
			sine = bestSine * 2.0 - sine;
		}

		return sine/bestSine;
	}


	// calculate dihedral sine
	/** param[heh_]: indexe of the shared edge
	    return dihedral sine
	**/
	double EdgeContract::dihedral_sine(VolumeMesh::TetraMesh::HalfEdgeHandle heh_)
	{
		VolumeMesh::TetraMesh::HalfFaceHandle hf1, hf2;
		hf1 = VolumeMesh::TetraMesh::HalfFaceHandle(heh_.idx()/3);
		hf2 = VolumeMesh::TetraMesh::HalfFaceHandle(mesh->mate_half_edge_handle(heh_).idx()/3);

		VolumeMesh::TetraMesh::Normal fn1, fn2;
		fn1 = mesh->normal(hf1);
		fn2 = mesh->normal(hf2);

		double sine, bestSine;
		bestSine = pow(3.0, 0.5) / 2.0;  // sin(60°)
		//sine = pow(1 - pow((n[fk] | n[fl]) / n[fk].norm() / n[fl].norm(), 2.0), 0.5);
		sine = pow(1 - pow((fn1 | fn2)/fn1.norm()/fn2.norm(),2.0),0.5);
		if ((sine / bestSine) > 1)
		{
			sine = bestSine * 2.0 - sine;
		}

		return sine/bestSine;
	}

	// jacobian value measure of a tetrahedron
	double EdgeContract::minimun_jacobian(VolumeMesh::Point v[4])
	{
		double jacobian = 1, temp;
		for (int i = 0; i < 4; i ++)
		{
			temp = fabs(jacobian_value(v,i));
			if (temp < jacobian)
				jacobian = temp;
		}
		return jacobian;
	}

	// a single corner's jacobian value
	double EdgeContract::jacobian_value(VolumeMesh::Point v_[4], int vIdx)
	{
		VolumeMesh::Point v[4];
		for (int i = 0; i < 4; i ++)
			v[i] = v_[(vIdx+i)%4];

		double jacobian;
		VolumeMesh::TetraMesh::Normal n[3], w[3], ideal[3];
		w[0] = VolumeMesh::TetraMesh::Normal(1.0, 0, 0);
		w[1] = VolumeMesh::TetraMesh::Normal(-0.5774, 1.1547, 0);
		w[2] = VolumeMesh::TetraMesh::Normal(-0.4082, -0.4082, 1.2247);

		n[0] = (v[1] - v[0]).normalize();
		n[1] = (v[2] - v[0]).normalize();
		n[2] = (v[3] - v[0]).normalize();

		for (int k = 0 ; k < 3; k ++)
		{
			for (int j = 0; j < 3; j ++)
			{
				ideal[j][k] = VolumeMesh::TetraMesh::Normal(n[0][k], n[1][k], n[2][k]) | w[j];
			}
		}
		ideal[0].normalize();
		ideal[1].normalize();
		ideal[2].normalize();

		jacobian = ideal[0] | (ideal[1] % ideal[2]);
		return jacobian;
	}

	// ratio of shortest side and longest side
	double EdgeContract::side_ratio(VolumeMesh::TetraHandle handle)
	{
		VolumeMesh::Point p[4];
		VolumeMesh::TetraMesh::HedronVertexIter hv_iter;
		double sides[6];
		double longSide = 0, shortSide = 0;
		hv_iter = mesh->hedron_vertex_iter(handle);
		p[0] = mesh->point(hv_iter.handle());
		++ hv_iter;
		p[1] = mesh->point(hv_iter.handle());
		++ hv_iter;
		p[2] = mesh->point(hv_iter.handle());
		++ hv_iter;
		p[3] = mesh->point(hv_iter.handle());

		sides[0] = (p[1] - p[0]).norm();
		sides[1] = (p[2] - p[0]).norm();
		sides[2] = (p[3] - p[0]).norm();
		sides[3] = (p[2] - p[1]).norm();
		sides[4] = (p[3] - p[2]).norm();
		sides[5] = (p[1] - p[3]).norm();

		longSide = shortSide = sides[0];
		for (int i = 1; i < 6; i++)
		{
			if (longSide < sides[i])
				longSide = sides[i];
			if (shortSide > sides[i])
				shortSide = sides[i];
		}

		return shortSide/longSide;
	}
	//---------------------------------end hedron quality calculate---------------------------------//

	//收缩边一周的四面体网格的雅克比值计算
	void EdgeContract::jac_val_around_edge(VolumeMesh::TetraMesh::HalfEdgeHandle heh, std::vector<double> &jacValVec)
	{
		int startTetraIdx;  // 起始四面体，在edgeStar中的索引
		double val, t_val;
		VolumeMesh::Point pStart, pEnd;  // 
		int pIdxStart, pIdxEnd;
		double jacStart, jacEnd;
		std::vector<VolumeMesh::TetraHandle> edgeStar;
		std::vector<VolumeMesh::Point> tetraPoints;  // 存放 edgeStar 中四面体的顶点坐标
		VolumeMesh::Point tPoints[4];
		VolumeMesh::TetraMesh::HedronVertexIter hv_it;
		mesh->edge_star(heh, edgeStar);
		tetraPoints.resize(edgeStar.size()*4);

		//获取volumelength质量值最小的网格的索引 & edgeStar中所有四面体的顶点坐标
		t_val = 1;
		for (unsigned int i = 0; i < edgeStar.size(); ++i)
		{
			hv_it = mesh->hedron_vertex_iter(edgeStar[i]);
			tetraPoints[i*4]   = mesh->point(hv_it.handle()); ++ hv_it;
			tetraPoints[i*4+1] = mesh->point(hv_it.handle()); ++ hv_it;
			tetraPoints[i*4+2] = mesh->point(hv_it.handle()); ++ hv_it;
			tetraPoints[i*4+3] = mesh->point(hv_it.handle());

			tPoints[0] = tetraPoints[i*4];
			tPoints[1] = tetraPoints[i*4+1];
			tPoints[2] = tetraPoints[i*4+2];
			tPoints[3] = tetraPoints[i*4+3];
			val = volume_length(tPoints);
			if (val < t_val)
			{
				t_val = val;
				startTetraIdx = i;
			}
		}

		if (startTetraIdx < 0 || startTetraIdx >= (int)edgeStar.size())
		{
			std::cerr << "Tetrahedron searching error!" << std::endl;
			return;
		}

		// 求起始tetra两顶点处的jacobian value
		pStart = mesh->point(mesh->from_vertex_handle(heh));
		pEnd = mesh->point(mesh->to_vertex_handle(heh));
		for (int i = 0; i < 4; i ++)
		{
			if (tetraPoints[startTetraIdx*4+i]==pStart)
				pIdxStart = i;
			else if (tetraPoints[startTetraIdx*4+i]==pEnd)
				pIdxEnd = i;
		}

		if (pIdxStart<0 || pIdxStart>=4 || pIdxEnd<0 || pIdxEnd>=4)
		{
			std::cerr << "Point searching error!" << std::endl;
			return;
		}
		tPoints[0] = tetraPoints[startTetraIdx*4];
		tPoints[1] = tetraPoints[startTetraIdx*4+1];
		tPoints[2] = tetraPoints[startTetraIdx*4+2];
		tPoints[3] = tetraPoints[startTetraIdx*4+3];
		jacStart = fabs(jacobian_value(tPoints, pIdxStart));
		jacEnd = fabs(jacobian_value(tPoints, pIdxEnd));
		if (jacEnd < jacStart)
		{
			double t_d;
			t_d = jacStart;
			jacStart = jacEnd;
			jacEnd = t_d;
			VolumeMesh::Point t_p;
			t_p = pStart;
			pStart = pEnd;
			pEnd = t_p;
		}
		jacValVec.push_back(jacStart);
		jacValVec.push_back(jacEnd);

		// 求其他四面体的雅克比值
		for (unsigned int tidx = 1; tidx < edgeStar.size(); tidx++)
		{
			int currentTetra = (startTetraIdx+tidx) % edgeStar.size();
			for (int i = 0; i < 4; i ++)
			{
				if (tetraPoints[currentTetra*4+i]==pStart)
					pIdxStart = i;
				else if (tetraPoints[currentTetra*4+i]==pEnd)
					pIdxEnd = i;
			}

			if (pIdxStart<0 || pIdxStart>=4 || pIdxEnd<0 || pIdxEnd>=4)
			{
				std::cerr << "Point searching error!" << std::endl;
				return;
			}
			tPoints[0] = tetraPoints[currentTetra*4];
			tPoints[1] = tetraPoints[currentTetra*4+1];
			tPoints[2] = tetraPoints[currentTetra*4+2];
			tPoints[3] = tetraPoints[currentTetra*4+3];
			jacStart = fabs(jacobian_value(tPoints, pIdxStart));
			jacEnd = fabs(jacobian_value(tPoints, pIdxEnd));

			jacValVec.push_back(jacStart);
			jacValVec.push_back(jacEnd);
		}
	}
}