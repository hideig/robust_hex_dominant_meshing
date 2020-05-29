#include <TopologyProcess/TopologyFlip.h>
#include <stdlib.h>
void TopologyFlip::tetraHedronFlip(int flipType_, int attributeType_, int qualityMeasureType_)
{
	if (isOutPut)
	{
		THRESHOLD = 10e-4;  // 设置阈值
		//duration = 0;
		QueryPerformanceFrequency(&li_start);
		PCFreq = double(li_start.QuadPart);
		timeCounter = 0;
		// initially set
		flipType = flipType_;
		attributeType = attributeType_;
		qualityMeasureType = qualityMeasureType_;

		FlipTetraPairNum = 0;
		FlipTimes = 0;
		FlipTimes_succeeded = 0;
		FlipTimes_failed = 0;
		totalTetraNum = 0;

		// create record files' names
		std::string s_flipType, s_qualityType;
		std::string fileName, fileName2, fileName3, fileName_tc;
		fileName = "F:\\testResult_new\\svm_flip";
		fileName2 = "F:\\testResult_new\\Q_flip";
		fileName3 = "F:\\testResult_new\\flip";
		fileName_tc = "F:\\flip_time_consumption";
		outfile_tc.open(fileName_tc.c_str(), std::fstream::app);
		if (!outfile_tc)
			return;

		if (flipType == 1)
		{
			fileName += "22_";
			fileName2 += "22_";
			fileName3 += "22_";
			s_flipType = "FLIP22";
			outfile_tc << "flip22_";
		}
		if (flipType == 2)
		{
			fileName += "23_";
			fileName2 += "23_";
			fileName3 += "23_";
			s_flipType = "FLIP23";
			outfile_tc << "flip23_";
		}
		if (flipType == 3)
		{
			fileName += "32_";
			fileName2 += "32_";
			fileName3 += "32_";
			s_flipType = "FLIP32";
			outfile_tc << "flip32_";
		}
		if (flipType == 4)
		{
			fileName += "44_";
			fileName2 += "44_";
			fileName3 += "44_";
			s_flipType = "FLIP44";
			outfile_tc << "flip44_";
		}

		fileName += entityName;
		fileName2 += entityName;
		outfile_tc << entityName;

		if (qualityMeasureType & 0x0001)
		{
			fileName += "_VolLen";
			fileName2 += "_VolLen";
			fileName3 += "_VolLen";
			s_qualityType = "VOLUMLENGTH";
			outfile_tc << "_VolLen";
		}
		if (qualityMeasureType & 0x0010)
		{
			fileName += "_MinSin";
			fileName2 += "_MinSin";
			fileName3 += "_MinSin";
			s_qualityType = "MINIMUMSINE";
			outfile_tc << "_MinSin";
		}
		if (qualityMeasureType & 0x0100)
		{
			fileName += "_JacVal";
			fileName2 += "_JacVal";
			fileName3 += "_JacVal";
			s_qualityType = "JACOBIANVALUE";
			outfile_tc << "_JacVal";
		}

		outfile.open(fileName.c_str());
		outfile2.open(fileName2.c_str());
		outfile3.open(fileName3.c_str(), ofstream::app);

		if (!outfile || !outfile2 || !outfile3)
		{
			std::cerr<<"error: unable to open input file" ;
			return;
		}

		// outfile2 data save
		VolumeMesh::TetraMesh::HedronIter t_iter;
		VolumeMesh::TetraMesh::HedronVertexIter hv_iter;
		VolumeMesh::Point tpoint[4];
		// save flip type and quality type
		outfile2 << s_flipType << " " << s_qualityType << "\n";
		// save total tetrahedron number
		outfile2 << tetraMesh->size_tetrahedron() << "\n";
		// save original tetrahedron quality
		for (t_iter = tetraMesh->hedrons_begin(); t_iter != tetraMesh->hedrons_end(); ++ t_iter)
		{
			hv_iter = tetraMesh->hedron_vertex_iter(t_iter);
			tpoint[0] = tetraMesh->point(hv_iter.handle());
			++ hv_iter;
			tpoint[1] = tetraMesh->point(hv_iter.handle());
			++ hv_iter;
			tpoint[2] = tetraMesh->point(hv_iter.handle());
			++ hv_iter;
			tpoint[3] = tetraMesh->point(hv_iter.handle());

			outfile2 << getQualityValue(tpoint) << "\n";
		}

		// do flip
		tetras_NeedToFlip.clear();
		oldTetraQDS();  // 原始网格的质量分布
		if (flipType == 1 || flipType == 2)
		{
			flip_22_23();
		}
		else if (flipType == 3)
		{
			flip_32();
		}
		else if (flipType == 4)
		{
			flip_44();
		}
		newTetraQDS();  // 新网格的质量分布

		// outfile3 data save
		totalTetraNum = tetraMesh->size_tetrahedron();
		outfile3 << entityName << "\t" << totalTetraNum << "\t" << FlipTetraPairNum << "\t" 
			<< FlipTimes_succeeded << "\t" << FlipTimes_failed << "\t" << (FlipTimes_succeeded * 1.0)/FlipTimes_failed;
		for (int i = 0; i < 5; i++)
		{
			outfile3 << "\t" << oldQD[i];
		}
		for (int i = 0; i < 5; i++)
		{
			outfile3 << "\t" << newQD[i];
		}
		outfile3 << "\n";

		//outfile_tc << "\t" << duration/CLOCKS_PER_SEC << "\t" << (duration/CLOCKS_PER_SEC/FlipTetraPairNum) << "\n";
		outfile_tc << "\t" << timeCounter/PCFreq << "\t" << (timeCounter/PCFreq/FlipTetraPairNum) << "\n";

		outfile.close();
		outfile2.close();
		outfile3.close();
		outfile_tc.close();
	}
	else
	{
		//duration = 0;
		QueryPerformanceFrequency(&li_start);
		PCFreq = double(li_start.QuadPart);
		timeCounter = 0;
		// initially set
		flipType = flipType_;
		attributeType = attributeType_;
		qualityMeasureType = qualityMeasureType_;

		FlipTetraPairNum = 0;
		FlipTimes = 0;
		FlipTimes_succeeded = 0;
		FlipTimes_failed = 0;
		totalTetraNum = 0;

		std::string fileName_tc;
		fileName_tc = "F:\\flip_time_consumption";
		outfile_tc.open(fileName_tc.c_str(), std::fstream::app);
		if (!outfile_tc)
			return;

		tetras_NeedToFlip.clear();
		if (flipType == 1)
		{
			outfile_tc << "flip22_";
		}
		if (flipType == 2)
		{
			outfile_tc << "flip23_";
		}
		if (flipType == 3)
		{
			outfile_tc << "flip32_";
		}
		if (flipType == 4)
		{
			outfile_tc << "flip44_";
		}

		outfile_tc << entityName;

		if (qualityMeasureType & 0x0001)
		{
			outfile_tc << "_VolLen";
		}
		if (qualityMeasureType & 0x0010)
		{
			outfile_tc << "_MinSin";
		}
		if (qualityMeasureType & 0x0100)
		{
			outfile_tc << "_JacVal";
		}

		if (flipType == 1 || flipType == 2)
		{
			flip_22_23();
		}
		else if (flipType == 3)
		{
			flip_32();
		}
		else if (flipType == 4)
		{
			flip_44();
		}

		//outfile_tc << "\t" << duration/CLOCKS_PER_SEC << "\t" << (duration/CLOCKS_PER_SEC/FlipTetraPairNum) << "\n";
		outfile_tc << "\t" << timeCounter/PCFreq << "\t" << (timeCounter/PCFreq/FlipTetraPairNum) << "\n";
		outfile_tc.close();
	}
}

void TopologyFlip::flip_22_23()
{
	VolumeMesh::TetraMesh::HedronIter th_iter;
	VolumeMesh::TetraMesh::HedronFaceIter hf_iter;
	VolumeMesh::TetraMesh::HalfFaceHandle oppoHFace;
	VolumeMesh::TetraMesh::HedronHandle oppoTetra;
	std::vector<VolumeMesh::TetraMesh::HedronHandle> tempTetra;
	VolumeMesh::TetraMesh::HedronHandle flipTetra[2];      //flipping tetras
	VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2];  //shared half face of two flipping tetras
	VolumeMesh::TetraMesh::HalfFaceHandle coplanarFace[2]; //coplanar faces

	flippedTetraComb.clear();
	for (th_iter = tetraMesh->hedrons_begin(); th_iter != tetraMesh->hedrons_end(); ++ th_iter)
	{
		flipTetra[0] = th_iter.handle();
		//iterate every half face to do 2-3 or 2-2 flip
		for (hf_iter = tetraMesh->hedron_face_iter(th_iter); hf_iter; ++ hf_iter)
		{
			if (tetraMesh->has_opposite_half_face(hf_iter.handle()))
			{
				sharedHFace[0] = hf_iter.handle();
				oppoHFace = tetraMesh->opposite_half_face_handle(hf_iter.handle());
				oppoTetra = tetraMesh->hedron_handle(oppoHFace);
				flipTetra[1] = oppoTetra;
				sharedHFace[1] = oppoHFace;

				//check if the combination of the tetras have been flipped
				std::string tetraComb;
				tempTetra.clear();
				tempTetra.push_back(flipTetra[0]);
				tempTetra.push_back(flipTetra[1]);
				if (is_tetraComb_flipped(tempTetra, tetraComb))  continue;
				flippedTetraComb.insert(tetraComb);

				//check if the combination of the two tetras is convex 
				//start = clock();
				QueryPerformanceCounter(&li_start);
				bool isConvex;
				isConvex = is_convex(sharedHFace[0]);
				QueryPerformanceCounter(&li_finish);
				timeCounter += double(li_finish.QuadPart - li_start.QuadPart);
				//flinish = clock();
				//duration += (double)(flinish - start);
				if (!isConvex)  continue;

				//do flip23 or flip22
				flip_23(flipTetra, sharedHFace);
			}
		}
	}
}

void TopologyFlip::flip_23(VolumeMesh::TetraMesh::HedronHandle flipTetra[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2])
{
	VolumeMesh::Point p[4];
	VolumeMesh::TetraMesh::Normal n[3];
	VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
	VolumeMesh::TetraMesh::HedronFaceIter hf_iter;
	VolumeMesh::TetraMesh::HalfFaceVertexIter fv_iter;
	VolumeMesh::TetraMesh::HedronVertexIter hv_iter;
	TetraHedronEntity * newTetra;
	std::vector<TetraHedronEntity> flippedTetra, originalTetra;
	double volume;
	int counter;

	//start = clock();
	QueryPerformanceCounter(&li_start);
	for (int i = 0; i < 2; i++)
	{
		hv_iter = tetraMesh->hedron_vertex_iter(flipTetra[i]);
		p[0] = tetraMesh->point(hv_iter.handle()); ++ hv_iter;
		p[1] = tetraMesh->point(hv_iter.handle()); ++ hv_iter;
		p[2] = tetraMesh->point(hv_iter.handle()); ++ hv_iter;
		p[3] = tetraMesh->point(hv_iter.handle());
		originalTetra.push_back(TetraHedronEntity(p));
	}

	//find the point in sharedHFace[1] but not in sharedHFace[0]
	fe_iter = tetraMesh->face_half_edge_iter(sharedHFace[1]);
	p[0] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle()))));

	//flip two tetra into three 
	counter = 0;
	newTetra = new TetraHedronEntity[3];
	for (hf_iter = tetraMesh->hedron_face_iter(flipTetra[0]); hf_iter; ++ hf_iter)
	{
		if (hf_iter.handle() == sharedHFace[0]) continue;
		fv_iter = tetraMesh->half_face_vertex_iter(hf_iter.handle());
		p[1] = tetraMesh->point(fv_iter.handle()); ++ fv_iter;
		p[2] = tetraMesh->point(fv_iter.handle()); ++ fv_iter;
		p[3] = tetraMesh->point(fv_iter.handle());
		newTetra[counter++].setPoint(p);
	}

	//calculate the volumes of new tetras
	for (int i = 0; i < 3; i ++)
	{
		newTetra[i].getPoint(p);
		n[0] = p[1] - p[0];
		n[1] = p[2] - p[0];
		n[2] = p[3] - p[0];
		volume = ((n[0] % n[1]) | n[2]) / 6;
		//delete the zero volume tetra(2-2 flip)
		if (fabs(volume) > 1e-6)
			flippedTetra.push_back(newTetra[i]);
	}

	if ((is_optimized(flippedTetra, originalTetra)) &&
		((flipType == 1 && flippedTetra.size() == 2) || (flipType == 2 && flippedTetra.size() == 3)))
	{
		tetras_NeedToFlip.insert(flipTetra[0]);
		tetras_NeedToFlip.insert(flipTetra[1]);

		// find two old tetras
		VolumeMesh::TetraMesh::HedronIter th_iter;
		for (th_iter = tetraMesh->hedrons_begin(); th_iter != tetraMesh->hedrons_end() && th_iter.handle() != flipTetra[0]; ++ th_iter);
		for (th_iter = tetraMesh->hedrons_begin(); th_iter != tetraMesh->hedrons_end() && th_iter.handle() != flipTetra[1]; ++ th_iter);
	}

	delete [] newTetra;
	QueryPerformanceCounter(&li_finish);
	timeCounter += double(li_finish.QuadPart - li_start.QuadPart);
	//flinish = clock();
	//duration += (double)(flinish - start);

	//save data
	if (isOutPut)
	{
		if (flipType == 1 && flippedTetra.size() == 2)
		{
			dataSave_flip22(originalTetra, flippedTetra, flipTetra, sharedHFace);
		}
		else if (flipType == 2 && flippedTetra.size() == 3)
		{
			dataSave_flip23(originalTetra, flippedTetra, flipTetra, sharedHFace);
		}
	}
}


void TopologyFlip::flip_32()
{
	VolumeMesh::TetraMesh::HedronIter h_iter;
	VolumeMesh::TetraMesh::HedronHalfEdgeIter hhe_iter;
	VolumeMesh::TetraMesh::HedronVertexIter hv_iter;
	std::vector<VolumeMesh::TetraMesh::HedronHandle> edgeStar;
	std::vector<TetraHedronEntity> oldtetra, flipTetra_32;
	std::set<VolumeMesh::Point> pSet;
	std::set<VolumeMesh::Point>::iterator pSet_it;
	VolumeMesh::Point tp, p[5], newTetraPoint1[4], newTetraPoint2[4];
	VolumeMesh::TetraMesh::Normal newTetraEdge1[3], newTetraEdge2[3];
	double volume1, volume2;

	//遍历所有的四面体，若有三个网格共享一条边且该边非边界边，进行3-2flip
	flippedTetraComb.clear();
	for (h_iter = tetraMesh->hedrons_begin(); h_iter != tetraMesh->hedrons_end(); ++ h_iter)
	{
		for (hhe_iter = tetraMesh->hedron_half_edge_iter(h_iter.handle()); hhe_iter; ++hhe_iter)
		{
			//start = clock();
			QueryPerformanceCounter(&li_start);
			edgeStar.clear();
			tetraMesh->edge_star(hhe_iter.handle(),edgeStar);
			QueryPerformanceCounter(&li_finish);
			timeCounter += double(li_finish.QuadPart - li_start.QuadPart);
			//flinish = clock();
			//duration += (double)(flinish - start); 

			if (edgeStar.size() == 3 && !tetraMesh->is_boundary(hhe_iter.handle()))
			{
				//start = clock();
				QueryPerformanceCounter(&li_start);
				//check if the combination of the tetras have been flipped
				std::string tetraComb;
				if (is_tetraComb_flipped(edgeStar, tetraComb))
				{
					continue;
				}
				flippedTetraComb.insert(tetraComb);

				//do 3-2 flip
				p[0] = tetraMesh->point(tetraMesh->from_vertex_handle(hhe_iter.handle()));
				p[1] = tetraMesh->point(tetraMesh->to_vertex_handle(hhe_iter.handle()));

				//get the other three points and old tetras
				pSet.clear();
				oldtetra.clear();
				for (unsigned int i = 0; i < edgeStar.size(); i ++)
				{
					VolumeMesh::Point tempPoint[4];
					int j = 0;
					for (hv_iter = tetraMesh->hedron_vertex_iter(edgeStar[i]); hv_iter; ++ hv_iter)
					{
						tp = tetraMesh->point(hv_iter.handle());
						if (tp != p[0] && tp != p[1])
						{
							pSet.insert(tp);
						}
						tempPoint[j++] = tp;
					}
					oldtetra.push_back(TetraHedronEntity(tempPoint));
				}

				if (pSet.size() == 3)
				{
					pSet_it = pSet.begin();
					p[2] = *pSet_it;
					++ pSet_it;
					p[3] = * pSet_it;
					++ pSet_it;
					p[4] = * pSet_it;
				}

				//get two new tetras points
				newTetraPoint1[0] = p[0];
				newTetraPoint1[1] = p[2];
				newTetraPoint1[2] = p[3];
				newTetraPoint1[3] = p[4];

				newTetraPoint2[0] = p[1];
				newTetraPoint2[1] = p[2];
				newTetraPoint2[2] = p[4];
				newTetraPoint2[3] = p[3];

				//check if there is the reversal tetra or zero volume tetra
				newTetraEdge1[0] = newTetraPoint1[1] - newTetraPoint1[0];
				newTetraEdge1[1] = newTetraPoint1[2] - newTetraPoint1[0];
				newTetraEdge1[2] = newTetraPoint1[3] - newTetraPoint1[0];

				newTetraEdge2[0] = newTetraPoint2[1] - newTetraPoint2[0];
				newTetraEdge2[1] = newTetraPoint2[2] - newTetraPoint2[0];
				newTetraEdge2[2] = newTetraPoint2[3] - newTetraPoint2[0];

				volume1 = (newTetraEdge1[0] % newTetraEdge1[1]) | newTetraEdge1[2];
				volume2 = (newTetraEdge2[0] % newTetraEdge2[1]) | newTetraEdge2[2];

				//get the new tetras, delete the zero volume tetra
				flipTetra_32.clear();
				if ((volume1 * volume2) == 0)
				{
					if (volume1 == 0)
					{
						if (volume2 < 0)
						{
							tp = newTetraPoint2[2];
							newTetraPoint2[2] = newTetraPoint2[3];
							newTetraPoint2[3] = tp;
						}
						flipTetra_32.push_back(TetraHedronEntity(newTetraPoint2));
					}
					else
					{
						if (volume1 < 0)
						{
							tp = newTetraPoint1[2];
							newTetraPoint1[2] = newTetraPoint1[3];
							newTetraPoint1[3] = tp;
						}
						flipTetra_32.push_back(TetraHedronEntity(newTetraPoint1));
					}
				}
				else if ((volume1 * volume2) > 0)
				{
					if (volume1 < 0 && volume2 < 0)
					{
						tp = newTetraPoint1[2];
						newTetraPoint1[2] = newTetraPoint1[3];
						newTetraPoint1[3] = tp;

						tp = newTetraPoint2[2];
						newTetraPoint2[2] = newTetraPoint2[3];
						newTetraPoint2[3] = tp;
					}
					flipTetra_32.push_back(TetraHedronEntity(newTetraPoint1));
					flipTetra_32.push_back(TetraHedronEntity(newTetraPoint2));
				}
				bool isOpt = is_optimized(flipTetra_32, oldtetra);  // timing needed
				QueryPerformanceCounter(&li_finish);
				timeCounter += double(li_finish.QuadPart - li_start.QuadPart);
				//flinish = clock();
				//duration += (double)(flinish - start); 

				if (!flipTetra_32.size())
					continue;

				// save data
				if (isOutPut)
				{
					dataSave_flip32(oldtetra, flipTetra_32, p, edgeStar);
				}

				//check if the quality is improved after 3-2 flip
				if (is_optimized(flipTetra_32, oldtetra))
				{
					for (unsigned int i = 0; i < edgeStar.size(); i++)
					{
						tetras_NeedToFlip.insert(edgeStar[i]);
					}
				}
			}
		}
	}
}

//there are two schemes of 4-4 flip
void TopologyFlip::flip_44()
{
	std::vector<TetraHedronEntity> oldTetra, newTetra;
	VolumeMesh::TetraMesh::HedronIter h_iter;
	VolumeMesh::TetraMesh::HedronHalfEdgeIter he_iter;
	std::vector<VolumeMesh::TetraMesh::HedronHandle> edgeStar;
	VolumeMesh::Point p[6];
	VolumeMesh::Point tetraPoint[4];
	std::vector<TetraHedronEntity> scheme1Tetra, scheme2Tetra, originalTetra;
	bool comp1, comp2;

	flippedTetraComb.clear();
	for (h_iter = tetraMesh->hedrons_begin(); h_iter != tetraMesh->hedrons_end(); ++ h_iter)
	{
		for (he_iter = tetraMesh->hedron_half_edge_iter(h_iter.handle()); he_iter; ++he_iter)
		{
			edgeStar.clear();
			tetraMesh->edge_star(he_iter.handle(), edgeStar);
			if (edgeStar.size() == 4 && !tetraMesh->is_boundary(he_iter.handle()))
			{
				//start = clock();
				QueryPerformanceCounter(&li_start);
				//check if the combination of the tetras have been flipped
				std::string tetraComb;
				if (is_tetraComb_flipped(edgeStar, tetraComb))
				{
					continue;
				}
				flippedTetraComb.insert(tetraComb);

				//do 4-4 flip
				p[0] = tetraMesh->point(tetraMesh->from_vertex_handle(he_iter.handle()));
				p[1] = tetraMesh->point(tetraMesh->to_vertex_handle(he_iter.handle()));
				p[2] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(he_iter.handle())));
				p[3] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(he_iter.handle()))));
				p[4] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(tetraMesh->radial_half_edge_handle(tetraMesh->mate_half_edge_handle(he_iter.handle()))))));
				p[5] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(tetraMesh->radial_half_edge_handle(he_iter.handle())))));

				//get originalTetra
				tetraPoint[0] = p[3];tetraPoint[1] = p[0];tetraPoint[2] = p[1];tetraPoint[3] = p[2];
				originalTetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[3];tetraPoint[1] = p[0];tetraPoint[2] = p[4];tetraPoint[3] = p[1];
				originalTetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[5];tetraPoint[1] = p[0];tetraPoint[2] = p[2];tetraPoint[3] = p[1];
				originalTetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[5];tetraPoint[1] = p[0];tetraPoint[2] = p[1];tetraPoint[3] = p[4];
				originalTetra.push_back(TetraHedronEntity(tetraPoint));

				//get scheme1Tetra
				tetraPoint[0] = p[0];tetraPoint[1] = p[2];tetraPoint[2] = p[4];tetraPoint[3] = p[3];
				scheme1Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[0];tetraPoint[1] = p[2];tetraPoint[2] = p[5];tetraPoint[3] = p[4];
				scheme1Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[1];tetraPoint[1] = p[2];tetraPoint[2] = p[3];tetraPoint[3] = p[4];
				scheme1Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[1];tetraPoint[1] = p[2];tetraPoint[2] = p[4];tetraPoint[3] = p[5];
				scheme1Tetra.push_back(TetraHedronEntity(tetraPoint));

				//get scheme2Tetra
				tetraPoint[0] = p[2];tetraPoint[1] = p[0];tetraPoint[2] = p[3];tetraPoint[3] = p[5];
				scheme2Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[2];tetraPoint[1] = p[1];tetraPoint[2] = p[5];tetraPoint[3] = p[3];
				scheme2Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[4];tetraPoint[1] = p[0];tetraPoint[2] = p[5];tetraPoint[3] = p[3];
				scheme2Tetra.push_back(TetraHedronEntity(tetraPoint));
				tetraPoint[0] = p[4];tetraPoint[1] = p[1];tetraPoint[2] = p[3];tetraPoint[3] = p[5];
				scheme2Tetra.push_back(TetraHedronEntity(tetraPoint));

				QueryPerformanceCounter(&li_finish);
				timeCounter += double(li_finish.QuadPart - li_start.QuadPart);
				//flinish = clock();
				//duration += (double)(flinish - start);

				//check if the quality is improved after flip and choose a scheme
				comp1 = is_optimized(scheme1Tetra, originalTetra);
				comp2 = is_optimized(scheme2Tetra, originalTetra);

				//if both scheme1 and scheme2 are not better than original, then not do flip
				if (!comp1 && !comp2)
				{
					continue;
				}

				for (unsigned int i = 0; i < edgeStar.size(); i++)
				{
					tetras_NeedToFlip.insert(edgeStar[i]);
				}
			}//end if (edgeStar is 4)
		}//end half_edge iteration
	}//end hedrons iteration 
}

//check if the quality of the tetra1 is better than tetra2
bool TopologyFlip::is_optimized(std::vector<TetraHedronEntity> tetra1, std::vector<TetraHedronEntity> tetra2)
{
	bool flag_vol_len = false, flag_min_sin = false, flag_jac_val = false;
	if (qualityMeasureType & 0x0001)
		flag_vol_len = true;
	if (qualityMeasureType & 0x0010)
		flag_min_sin = true;
	if (qualityMeasureType & 0x0100)
		flag_jac_val = true;

	double Q_Vol_Len_1 = 1, Q_Vol_Len_2 = 1;
	double Q_Min_Sin_1 = 1, Q_Min_Sin_2 = 1;
	double Q_Jac_Val_1 = 1, Q_Jac_Val_2 = 1;
	double quality_temp;
	VolumeMesh::Point p[4];
	for (unsigned int i = 0; i < tetra1.size(); i++)
	{
		tetra1[i].getPoint(p);

		quality_temp = volume_length(p);
		if (quality_temp < Q_Vol_Len_1)
			Q_Vol_Len_1 = quality_temp;

		quality_temp = minimum_sine(p);
		if (quality_temp < Q_Min_Sin_1)
			Q_Min_Sin_1 = quality_temp;

		quality_temp = minimun_jacobian(p);
		if (quality_temp < Q_Jac_Val_1)
			Q_Jac_Val_1 = quality_temp;
	}
	for (unsigned int i = 0; i < tetra2.size(); i ++)
	{
		tetra2[i].getPoint(p);

		quality_temp = volume_length(p);
		if (quality_temp < Q_Vol_Len_2)
			Q_Vol_Len_2 = quality_temp;

		quality_temp = minimum_sine(p);
		if (quality_temp < Q_Min_Sin_2)
			Q_Min_Sin_2 = quality_temp;

		quality_temp = minimun_jacobian(p);
		if (quality_temp < Q_Jac_Val_2)
			Q_Jac_Val_2 = quality_temp;
	}

	bool isOpt = false;
	if (flag_vol_len && (Q_Vol_Len_1 - Q_Vol_Len_2 > THRESHOLD )) isOpt = true;
	if (flag_min_sin && (Q_Min_Sin_1 - Q_Min_Sin_2 > THRESHOLD )) isOpt = true;
	if (flag_jac_val && (Q_Jac_Val_1 - Q_Jac_Val_2 > THRESHOLD )) isOpt = true;
	return isOpt;
}

//Volume-length measure, the best quality value is 1
double TopologyFlip::volume_length(VolumeMesh::Point v[4])
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
double TopologyFlip::minimum_sine(VolumeMesh::Point v[4])
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

// calculate dihedral sine
/** param[v]: points of tetrahedron
    param[k, l]: indexes of vertex of the shared edge
	**/
double TopologyFlip::dihedral_sine(VolumeMesh::Point v[4], int k, int l)
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

// jacobian value measure of a tetrahedron
double TopologyFlip::minimun_jacobian(VolumeMesh::Point v[4])
{
	double jacobian = 1, temp;
	VolumeMesh::Point cornerJv[4];

	for (int i = 0; i < 4; i ++)
	{
		cornerJv[0] = v[i];
		cornerJv[1] = v[(i+1)%4];
		cornerJv[2] = v[(i+2)%4];
		cornerJv[3] = v[(i+3)%4];
		temp = fabs(jacobian_value(cornerJv));

		if (temp < jacobian)
		{
			jacobian = temp;
		}
	}
	return jacobian;
}

// a single corner's jacobian value
double TopologyFlip::jacobian_value(VolumeMesh::Point v[4])
{
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

	jacobian = (ideal[0] % ideal[1]) | ideal[2];
	return jacobian;
}

// ||p1 - p2|| / |area|
double TopologyFlip::edge_area(VolumeMesh::Point len_p1, VolumeMesh::Point len_p2, std::vector<VolumeMesh::Point> polygon)
{
	double length;
	length =  (len_p2 - len_p1).norm();

	double area = 0;
	VolumeMesh::TetraMesh::Normal edge1, edge2;
	for (unsigned int i = 2; i < polygon.size(); i ++)
	{
		edge1 = polygon[i - 1] - polygon[0];
		edge2 = polygon[i] - polygon[0];

		area += pow((edge1 | edge1) * (edge2 | edge2) - pow(edge1 | edge2, 2), 0.5);
	}
	area /= 2.0;

	return length/area;
}


// |area| / (a*a + b*b + c*c)
double TopologyFlip::area_length(VolumeMesh::Point triangle[3])
{
	double area;
	double length = 0;
	VolumeMesh::TetraMesh::Normal n[3];
	n[0] = triangle[0] - triangle[2];
	n[1] = triangle[1] - triangle[2];
	n[2] = triangle[0] - triangle[1];

	area = (n[0] % n[1]).norm();
	for (int i = 0; i < 3; i ++)
	{
		for (int j = 0; j < 3; j ++)
		{
			length += n[i][j] * n[i][j];
		}
	}

	return 2 * pow(3, 0.5) * area / length;
}

/**if the combination of the tetras have been flipped return true else return false
   tetraComb can return the unique string of the tetras' combination
*/
bool TopologyFlip::is_tetraComb_flipped(std::vector<VolumeMesh::TetraMesh::HedronHandle> tetra, std::string &tetraComb)
{
	char c[8];
	std::string * s;
	std::set<std::string>::iterator it;
	s = new std::string[tetra.size()];
	tetraComb = "";
	sort(tetra.begin(), tetra.end());

	itoa(tetra[0].idx(), c, 10);
	s[0] = (std::string)c;
	tetraComb += s[0];
	for (unsigned int i = 1; i < tetra.size(); i++)
	{
		itoa(tetra[i].idx(), c, 10);
		s[i] = (std::string)c;
		tetraComb = tetraComb + " " + s[i];
	}

	it = find(flippedTetraComb.begin(), flippedTetraComb.end(), tetraComb);
	if (it == flippedTetraComb.end())
	{
		return false;
	}
	return true;
}

// check if the combination of two tetrahedron is convex(return true)
bool TopologyFlip::is_convex(VolumeMesh::TetraMesh::HalfFaceHandle oneOfsharedHalfFace)
{
	VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
	VolumeMesh::Point v[4];
	VolumeMesh::TetraMesh::Normal n[2];
	bool isConvex = true;

	for (fe_iter = tetraMesh->face_half_edge_iter(oneOfsharedHalfFace); fe_iter; ++fe_iter)
	{
		v[0] = tetraMesh->point(tetraMesh->from_vertex_handle(fe_iter.handle()));
		v[1] = tetraMesh->point(tetraMesh->to_vertex_handle(fe_iter.handle()));
		v[2] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle()))));
		v[3] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(tetraMesh->radial_half_edge_handle(fe_iter.handle())))));

		n[0] = (v[1] - v[0]) % (v[2] - v[0]);
		n[0].normalize();
		n[1] = (v[3] - v[0]) % (v[1] - v[0]);
		n[1].normalize();

		v[0] = tetraMesh->point(tetraMesh->from_vertex_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle())));
		v[1] =tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle())));
		if (((n[0] % n[1]) | (v[1] - v[0])) < 0)
		{
			isConvex = false;
		}
	}

	return isConvex;
}

void TopologyFlip::dataSave_flip22(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra, 
								   VolumeMesh::TetraMesh::HedronHandle flipTetraH[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2])
{
	double vl_1, vl_2;             // volume-length values of two originalTetras
	double areaLen1, areaLen2;     // area length values of two coplanar faces
	double ratio;                  // ratio of two flip edges
	double vertexJcb1, vertexJcb2; // jacobian values of two top vertex
	double tetraJcb1, tetraJcb2;   // jacobian values of two tetrahedrons
	double sin1[3], sin2[3];       // three dihedral sins around the two top vertex
	VolumeMesh::Point p[5], triangle[3], tempv[4];
	VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
	VolumeMesh::TetraMesh::HalfFaceVertexIter fv_iter;
	std::vector<VolumeMesh::Point> polygon;
	int attrIdx = 0;

	// data save
	for (fe_iter = tetraMesh->face_half_edge_iter(sharedHFace[0]); fe_iter; ++fe_iter)
	{
		// get all five points of two tetrahedrons
		p[0] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle()))));
		p[1] = tetraMesh->point(tetraMesh->from_vertex_handle(fe_iter.handle()));
		p[2] = tetraMesh->point(tetraMesh->to_vertex_handle(fe_iter.handle()));
		p[3] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(fe_iter.handle())));
		p[4] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(tetraMesh->radial_half_edge_handle(fe_iter.handle())))));

		VolumeMesh::TetraMesh::Normal n[2];
		n[0] = (p[0] - p[1]) % (p[2] - p[1]);
		n[0].normalize();
		n[1] = (p[1] - p[4]) % (p[2] - p[4]);
		n[1].normalize();

		// find the coplanar faces
		if ((n[0] - n[1]).norm() < 1e-4)
		{
			// output label (classification)
			if (is_optimized(flippedTetra, originalTetra))
			{
				outfile << "1 ";
				outfile2 << "1 ";
				++ FlipTimes_succeeded;

				std::vector< VolumeMesh::TetraMesh::HedronHandle> flipTetraHV;
				flipTetraHV.push_back(flipTetraH[0]);
				flipTetraHV.push_back(flipTetraH[1]);
				updateTetraQuality(originalTetra, flippedTetra,flipTetraHV);   // update tetra quality
			}
			else
			{
				outfile << "0 ";
				outfile2 << "0 ";
				++ FlipTimes_failed;
				++ FlipTimes;
			}

			// flip data save
			++ FlipTimes;
			++ FlipTetraPairNum;
			outfile2 << flipTetraH[0].idx() << " " << flipTetraH[1].idx();
			// flipped tetras' quality
			for (unsigned int i = 0; i < flippedTetra.size(); i++)
			{
				flippedTetra[i].getPoint(tempv);
				outfile2 << " " << getQualityValue(tempv);
			}
			outfile2 << "\n";

			// svm data save
			/*output attributes*/
			// attribute volume_length
			if (attributeType & 0x0001)
			{
				originalTetra[0].getPoint(tempv);
				vl_1 = volume_length(tempv);
				originalTetra[1].getPoint(tempv);
				vl_2 = volume_length(tempv);

				if (vl_1 < vl_2)
				{
					double temp = vl_1;
					vl_1 = vl_2;
					vl_2 = temp;

					VolumeMesh::Point temp_point;
					temp_point = p[4];
					p[4] = p[0];
					p[0] = temp_point;

					temp_point = p[1];
					p[1] = p[2];
					p[2] = temp_point;
				}
				outfile << (++attrIdx) << ":" << vl_1 << " ";
				outfile << (++attrIdx) << ":" << vl_2 << " ";
			}

			// attribute jac_corner
			if (attributeType & 0x0010)
			{
				// the points are ordered
				tempv[0] = p[0]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
				vertexJcb1 = fabs(jacobian_value(tempv));
				tempv[0] = p[4]; tempv[1] = p[1]; tempv[2] = p[3]; tempv[3] = p[2];
				vertexJcb2 = fabs(jacobian_value(tempv));

				if (vertexJcb1 < vertexJcb2)
				{
					double tempJcb = vertexJcb1;
					vertexJcb1 = vertexJcb2;
					vertexJcb2 = tempJcb;

					VolumeMesh::Point temp_point;
					temp_point = p[4];
					p[4] = p[0];
					p[0] = temp_point;
				}
				outfile << (++attrIdx) << ":" << vertexJcb1 << " ";
				outfile << (++attrIdx) << ":" << vertexJcb2 << " ";
			}

			// attribute dihedral sine
			if (attributeType & 0x0040)
			{
				// the points are ordered
				tempv[0] = p[0]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
				for (int i = 0; i < 3; i ++)
				{
					sin1[i] = dihedral_sine(tempv, 0, i+1);
				}
				tempv[0] = p[4]; tempv[1] = p[1]; tempv[2] = p[3]; tempv[3] = p[2];
				for (int i = 0; i < 3; i ++)
				{
					sin2[i] = dihedral_sine(tempv, 0, i+1);
				}

				// sort
				for (int i = 0; i < 2; i++)
				{
					for (int j = i + 1; j < 3; j ++)
					{
						if (sin1[i] < sin1[j])
						{
							double tempVL = sin1[i];
							sin1[i] = sin1[j];
							sin1[j] = tempVL;

							tempVL = sin2[i];
							sin2[i] = sin2[j];
							sin2[j] = tempVL;
						}
					}
				}

				// output data
				for (int i = 0; i < 3; i ++)
					outfile << (++attrIdx) << ":" << sin1[i] << " ";
				for (int i = 0; i < 3; i ++)
					outfile << (++attrIdx) << ":" << sin2[i] << " ";
			}

			// attribute jac_tetra
			if (attributeType & 0x0020)
			{
				tempv[0] = p[0]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
				tetraJcb1 = fabs(minimun_jacobian(tempv));
				tempv[0] = p[4]; tempv[1] = p[1]; tempv[2] = p[3]; tempv[3] = p[2];
				tetraJcb2 = fabs(minimun_jacobian(tempv));
				outfile << (++attrIdx) << ":" << tetraJcb1 << " ";
				outfile << (++attrIdx) << ":" << tetraJcb2 << " ";
			}

			// attribute edge_ratio
			if (attributeType & 0x0008)
			{
				if ((p[0] - p[4]).norm() != 0)
					ratio = (p[2] - p[1]).norm() / (p[0] - p[4]).norm();
				else
					ratio = 0;
				outfile << (++attrIdx) << ":" << ratio << " ";
			}

			// attribute area_length
			if (attributeType & 0x0002)
			{
				triangle[0] = p[1];
				triangle[1] = p[0];
				triangle[2] = p[2];
				areaLen1 = area_length(triangle);
				triangle[0] = p[1];
				triangle[1] = p[2];
				triangle[2] = p[4];
				areaLen2 = area_length(triangle);
				outfile << (++attrIdx) << ":" << areaLen1 << " ";
				outfile << (++attrIdx) << ":" << areaLen2 << " ";
			}
			outfile << "\n";
			break;
		}
	}
}

void TopologyFlip::dataSave_flip23(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra, 
								   VolumeMesh::TetraMesh::HedronHandle flipTetraH[2], VolumeMesh::TetraMesh::HalfFaceHandle sharedHFace[2])
{
	double vl_1, vl_2;              // volume-length values of two originalTetras
	double vertexJcb1, vertexJcb2;  // jacobian values of two top vertex
	double tetraJcb1, tetraJcb2;    // jacobian values of two tetrahedrons
	double edgeArea;                // edge-area value
	VolumeMesh::Point p[5], tempv[4], triangle[3];
	VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
	VolumeMesh::TetraMesh::HalfFaceVertexIter fv_iter;
	std::vector<VolumeMesh::Point> polygon;
	int attrIdx = 0;

	// output label (classification)
	if (is_optimized(flippedTetra, originalTetra))
	{
		outfile << "1 ";
		outfile2 << "1 ";
		++ FlipTimes_succeeded;

		std::vector< VolumeMesh::TetraMesh::HedronHandle> flipTetraHV;
		flipTetraHV.push_back(flipTetraH[0]);
		flipTetraHV.push_back(flipTetraH[1]);
		updateTetraQuality(originalTetra, flippedTetra,flipTetraHV);   // update tetra quality
	}
	else
	{
		outfile << "0 ";
		outfile2 << "0 ";
		++ FlipTimes_failed;
		++ FlipTimes;
	}

	// flip data save
	++ FlipTimes;
	++ FlipTetraPairNum;
	outfile2 << flipTetraH[0].idx() << " " << flipTetraH[1].idx();
	// flipped tetras' quality
	for (unsigned int i = 0; i < flippedTetra.size(); i++)
	{
		flippedTetra[i].getPoint(tempv);
		outfile2 << " " << getQualityValue(tempv);
	}
	outfile2 << "\n";

	// svm data save
	// get all five points of two tetrahedrons
	fe_iter = tetraMesh->face_half_edge_iter(sharedHFace[0]);
	p[0] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(fe_iter.handle()))));
	p[1] = tetraMesh->point(tetraMesh->from_vertex_handle(fe_iter.handle()));
	p[2] = tetraMesh->point(tetraMesh->to_vertex_handle(fe_iter.handle()));
	p[3] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(fe_iter.handle())));
	p[4] = tetraMesh->point(tetraMesh->to_vertex_handle(tetraMesh->next_half_edge_handle(tetraMesh->mate_half_edge_handle(tetraMesh->radial_half_edge_handle(fe_iter.handle())))));

	/*output attributes*/
	// attribute volume_length
	if (attributeType & 0x0001)
	{
		originalTetra[0].getPoint(tempv);
		vl_1 = volume_length(tempv);
		originalTetra[1].getPoint(tempv);
		vl_2 = volume_length(tempv);

		if (vl_1 < vl_2)
		{
			double temp = vl_1;
			vl_1 = vl_2;
			vl_2 = temp;

			VolumeMesh::Point temp_point;
			temp_point = p[4];
			p[4] = p[0];
			p[0] = temp_point;
		}
		outfile << (++attrIdx) << ":" << vl_1 << " ";
		outfile << (++attrIdx) << ":" << vl_2 << " ";
	}

	// attribute jac_corner
	if (attributeType & 0x0010)
	{
		tempv[0] = p[0]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
		vertexJcb1 = fabs(jacobian_value(tempv));
		tempv[0] = p[4]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
		vertexJcb2 = fabs(jacobian_value(tempv));

		if (vertexJcb1 < vertexJcb2)
		{
			double tempJcb = vertexJcb1;
			vertexJcb1 = vertexJcb2;
			vertexJcb2 = tempJcb;

			VolumeMesh::Point temp_point;
			temp_point = p[4];
			p[4] = p[0];
			p[0] = temp_point;
		}
		outfile << (++attrIdx) << ":" << vertexJcb1 << " ";
		outfile << (++attrIdx) << ":" << vertexJcb2 << " ";
	}

	// attribute jac_tetra
	if (attributeType & 0x0020)
	{
		tempv[0] = p[0]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
		tetraJcb1 = fabs(minimun_jacobian(tempv));
		tempv[0] = p[4]; tempv[1] = p[1]; tempv[2] = p[2]; tempv[3] = p[3];
		tetraJcb2 = fabs(minimun_jacobian(tempv));
		outfile << (++attrIdx) << ":" << tetraJcb1 << " ";
		outfile << (++attrIdx) << ":" << tetraJcb2 << " ";
	}

	// attribute edge-area-ratio
	if (attributeType & 0x0004)
	{
		polygon.push_back(p[1]);
		polygon.push_back(p[2]);
		polygon.push_back(p[3]);
		edgeArea = edge_area(p[0], p[4], polygon);
		outfile << (++attrIdx) << ":" << edgeArea << " ";
	}
	outfile << "\n";
}

void TopologyFlip::dataSave_flip32(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra, 
								   VolumeMesh::Point p[5], std::vector<VolumeMesh::TetraMesh::HedronHandle> edgeStar)
{
	double vl[3]; // volume-length values of three originalTetras
	double edgeArea;   // edge-area value
	double areaLen[3]; // area length values of two coplanar faces
	VolumeMesh::Point tempv[4], triangle[3];
	VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
	VolumeMesh::TetraMesh::HalfFaceVertexIter fv_iter;
	std::vector<VolumeMesh::Point> polygon;
	int attrIdx = 0;

	// output label (classification)
	if (is_optimized(flippedTetra, originalTetra))
	{
		outfile << "1 ";
		outfile2 << "1";
		++ FlipTimes_succeeded;
		updateTetraQuality(originalTetra, flippedTetra,edgeStar);   // update tetra quality
	}
	else
	{
		outfile << "0 ";
		outfile2 << "0";
		++ FlipTimes_failed;
		++ FlipTimes;
	}

	// flip data save
	++ FlipTimes;
	++ FlipTetraPairNum;
	// idx
	for (unsigned int i = 0; i < edgeStar.size() ; i++)
	{
		outfile2 << " " << edgeStar[i].idx();
	}
	// flipped tetras' quality
	for (unsigned int i = 0; i < flippedTetra.size(); i++)
	{
		flippedTetra[i].getPoint(tempv);
		outfile2 << " " << getQualityValue(tempv);
	}
	outfile2 << "\n";

	// svm data save
	// attribute volume-length values
	if (attributeType & 0x0001)
	{
		originalTetra[0].getPoint(tempv);
		vl[0] = volume_length(tempv);
		originalTetra[1].getPoint(tempv);
		vl[1] = volume_length(tempv);
		originalTetra[2].getPoint(tempv);
		vl[2] = volume_length(tempv);

		// sort the volume-length values
		for (int i = 0; i < 2; i++)
		{
			for (int j = i + 1; j < 3; j ++)
			{
				if (vl[i] < vl[j])
				{
					double tempVL = vl[i];
					vl[i] = vl[j];
					vl[j] = tempVL;

					TetraHedronEntity tempTetra = originalTetra[i];
					originalTetra[i] = originalTetra[j];
					originalTetra[j] = tempTetra;
				}
			}
		}
		outfile << (++attrIdx) << ":" << vl[0] << " ";
		outfile << (++attrIdx) << ":" << vl[1] << " ";
		outfile << (++attrIdx) << ":" << vl[2] << " ";
	}

	// attribute edge-area value
	if (attributeType & 0x0004)
	{
		for (int i = 2; i < 5; i ++)
		{
			polygon.push_back(p[i]);
		}
		edgeArea = edge_area(p[0], p[1], polygon);
		outfile << (++attrIdx) << ":" << edgeArea << " ";
	}

	// attribute area-length value
	if (attributeType & 0x0002)
	{
		for (int i = 0; i < 3; i ++)
		{
			std::set<VolumeMesh::Point> pointSet;
			std::set<VolumeMesh::Point>::iterator pointSet_it;
			originalTetra[i].getPoint(tempv);
			for (int j = 0; j < 4; j ++)
			{
				pointSet.insert(tempv[j]);
			}
			originalTetra[(i+1)%3].getPoint(tempv);

			int count = 0;
			for (int j = 0; j < 4; j ++)
			{
				if (count > 3)
					return;

				pointSet_it = find(pointSet.begin(), pointSet.end(), tempv[j]);
				if (pointSet_it != pointSet.end())
				{
					triangle[count ++] = tempv[j];
				}
			}
			areaLen[i] = area_length(triangle);
		}
		outfile << (++attrIdx) << ":" << areaLen[0] << " ";
		outfile << (++attrIdx) << ":" << areaLen[1] << " ";
		outfile << (++attrIdx) << ":" << areaLen[2] << " ";
	}
	outfile <<"\n";
}

void TopologyFlip::oldTetraQDS()
{
	tetraQuality.clear();
	tetraQuality.resize(tetraMesh->size_tetrahedron());
	oldQD[0] = oldQD[1] = oldQD[2] = oldQD[3] = oldQD[4] = 0;

	VolumeMesh::TetraMesh::HedronIter h_it;
	VolumeMesh::TetraMesh::HedronVertexIter hv_it;
	VolumeMesh::Point point[4];
	double qVal;
	// calculate every tetrahedron's quality value
	for (h_it = tetraMesh->hedrons_begin(); h_it != tetraMesh->hedrons_end(); ++ h_it)
	{
		hv_it = tetraMesh->hedron_vertex_iter(h_it);
		point[0] = tetraMesh->point(hv_it.handle());
		++ hv_it;
		point[1] = tetraMesh->point(hv_it.handle());
		++ hv_it;
		point[2] = tetraMesh->point(hv_it.handle());
		++ hv_it;
		point[3] = tetraMesh->point(hv_it.handle());

		// get quality value of the tetrahedron
		qVal = getQualityValue(point);
		tetraQuality[h_it.handle().idx()] = qVal;

		// quality distribute
		QualityDistrbtStatistic(qVal, oldQD);
	}
}

void TopologyFlip::newTetraQDS()
{
	newQD[0] = newQD[1] = newQD[2] = newQD[3] = newQD[4] = 0;
	for (unsigned int i = 0; i < tetraQuality.size(); i++)
	{
		QualityDistrbtStatistic(tetraQuality[i], newQD);
	}
}

void TopologyFlip::updateTetraQuality(std::vector<TetraHedronEntity> originalTetra, std::vector<TetraHedronEntity> flippedTetra,
						std::vector<VolumeMesh::TetraMesh::HedronHandle> flipTetraHV)
{
	double *qVal_old = new double[originalTetra.size()];
	double *qVal_new = new double[flippedTetra.size()];

	// get original and flipped tetrahedrons' quality values
	VolumeMesh::Point p[4];
	for (unsigned int i = 0; i < originalTetra.size(); i++)
	{
		originalTetra[i].getPoint(p);
		qVal_old[i] = getQualityValue(p);
	}
	for (unsigned int i = 0; i < flippedTetra.size(); i++)
	{
		flippedTetra[i].getPoint(p);
		qVal_new[i] = getQualityValue(p);
	}

	// sort
	for (int i = 1; i < originalTetra.size(); i++)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			if (qVal_old[i] < qVal_old[j])
			{
				double td = qVal_old[i];
				qVal_old[i] = qVal_old[j];
				qVal_old[j] = td;

				VolumeMesh::TetraMesh::HedronHandle thh = flipTetraHV[i];
				flipTetraHV[i] = flipTetraHV[j];
				flipTetraHV[j] = thh;

			}
		}
	}
	sort(&qVal_new[0], &qVal_new[flippedTetra.size()]);

	// update
	if (flipType == 1) // flip22
	{
		tetraQuality[flipTetraHV[0].idx()] = qVal_new[0];
		tetraQuality[flipTetraHV[1].idx()] = qVal_new[1];
	}
	else if (flipType == 2) // flip23
	{
		tetraQuality[flipTetraHV[0].idx()] = qVal_new[0];
		tetraQuality[flipTetraHV[1].idx()] = qVal_new[2];
		tetraQuality.push_back(qVal_new[1]);
	}
	else if (flipType == 3) // flip32
	{
		tetraQuality[flipTetraHV[0].idx()] = -1;
		tetraQuality[flipTetraHV[1].idx()] = -1;
		tetraQuality[flipTetraHV[2].idx()] = -1;
		tetraQuality.push_back(qVal_new[0]);
		tetraQuality.push_back(qVal_new[1]);
	}
}

double TopologyFlip::getQualityValue(VolumeMesh::Point p[4])
{
	double qVal;
	if (qualityMeasureType & 0x0001)
	{
		qVal = volume_length(p);
	}
	if (qualityMeasureType & 0x0010)
	{
		qVal = minimum_sine(p);
	}
	if (qualityMeasureType & 0x0100)
	{
		qVal = minimun_jacobian(p);
	}
	return qVal;
}

void TopologyFlip::QualityDistrbtStatistic(double qVal, int * arrayQD)
{
	if (qVal >= 0 && qVal < 0.1)
	{
		++ arrayQD[0];
	}
	else if (qVal >= 0.1 && qVal < 0.3)
	{
		++ arrayQD[1];
	}
	else if (qVal >= 0.3 && qVal < 0.5)
	{
		++ arrayQD[2];
	}
	else if (qVal >= 0.5 && qVal < 0.7)
	{
		++ arrayQD[3];
	}
	else if (qVal >= 0.7 && qVal <= 1.0)
	{
		++ arrayQD[4];
	}
}