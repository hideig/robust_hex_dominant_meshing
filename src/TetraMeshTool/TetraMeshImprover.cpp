#include <TetraMeshTool/TetraMeshImprover.h>
namespace VolumeMesh
{
	/* initialize the improver */
	void TetraMeshImprover::initialize()
	{		
		/*control parameters initialization*/

		qualitymetric = QUAL_MINSINE;

		edgeremoval = 1;
		singlefaceremoval = 1;
		multifaceremoval = 1;

		/*tools initialization*/
		vsmoother.setMesh(mesh);
		vnonsmoother.setTmesh(mesh);
		vinserter.setMesh(mesh);
		quadric.setMesh(mesh);
		toperator.setMesh(mesh);
		qualcalculator.setMesh(mesh);

		vsmoother.setQualityMetric(qualitymetric);
		vinserter.setQualityMetric(qualitymetric);
		quadric.setQualityMetric(qualitymetric);
		toperator.setQualityMetric(qualitymetric);
		qualcalculator.setQualityMetric(qualitymetric);

		vsmoother.setStats(&stats);
		toperator.setStats(&stats);
		vsmoother.setImproveBehavior(&improvebehave);
		toperator.setImproveBehavior(&improvebehave);

		vsmoother.setJournal(&journals);
		vsmoother.setQuadricContainer(&quadric);
		vsmoother.setQualityCalculator(&qualcalculator);
		
		toperator.setJournal(&journals);
		toperator.setQualityCalculator(&qualcalculator);


		/*tools' control parameter initialization*/

		toperator.edgeremoval = edgeremoval;
		toperator.singlefaceremoval = singlefaceremoval;
		toperator.multifaceremoval = multifaceremoval;

		parseimprovecommandline(&improvebehave);
		parseimprovestateline(&stats);
	}

	void TetraMeshImprover::parseimprovecommandline(ImproveBehavior *b)
	{
		/* set default values for all options */
		b->qualmeasure = QUAL_MINSINE;
		/*
		b->qualmeasure = QUALWARPEDMINSINE;
		b->qualmeasure = QUALRADIUSRATIO;
		b->qualmeasure = QUALVLRMS3RATIO;
		*/
		b->sinewarpfactor = 0.75;
		b->usequadrics = 1;
		b->quadricoffset = 0.8;
		b->quadricscale = 300.0;

		b->nonsmooth = 1;
		b->facetsmooth = 1;
		b->segmentsmooth = 1;
		b->fixedsmooth = 0;

		b->edgeremoval = 1;
		b->boundedgeremoval = 0;
		b->singlefaceremoval = 1;
		b->multifaceremoval = 1;
		b->flip22 = 1;
		b->jflips = 1;

		b->insertthreshold = 0.035;
		b->enableinsert = 1;
		b->insertfacet = 1;
		b->insertbody = 1;
		b->insertsegment = 1;
		b->cavityconsiderdeleted = 0;
		b->cavdepthlimit = 6;

		/* Special Insertion options */
		b->enablespecialinsert = 1;
		b->specialinsertthreshold = 0.3;
		b->specialinsertbody = 1;
		b->specialinsertfacet = 1;

		b->edgecontraction = 1;

		b->minstepimprovement = 1.0e-4;
		b->mininsertionimprovement = 1.0e-3;
		b->maxinsertquality[0] = SINE40;
		b->maxinsertquality[1] = 0.7;
		b->maxinsertquality[2] = 0.7; 
		b->maxinsertquality[5] = SINE40;

		b->anisotropic = 0;
		b->tensor = 0;
		b->tensorb = 0;
		b->tensorblend = 1.0;

		b->sizing = 0;
		b->sizingpass = 0;
		b->targetedgelength = 0.0;
		b->longerfactor = 2.0;
		b->shorterfactor = 0.5;

		//b->verbosity = 1;
		b->usecolor = 0;

		b->outputandquit = 0;

		b->minsineout = 1;
		b->minangout = 0;
		b->maxangout = 0;
		b->vlrmsout = 0;
		b->nrrout = 0;

		strcpy(b->fileprefix, "");
		/*
		strcpy(b->fileprefix, "default");
		*/

		b->animate = 0;
		b->timeseries = 0;

		b->goalanglemin = 90.0;
		b->goalanglemax = 90.0;

		/* loop through each argument */
		//for (i=0; i<argc; i++)
		//{
		//	/* check if this is L */
		//	if (strcmp(argv[i], "-L") == 0)
		//	{
		//		/* set verbosity */
		//		improvebehave.verbosity = atoi(argv[i+1]);
		//	}

		//	/* check if this is F */
		//	if (strcmp(argv[i], "-F") == 0)
		//	{
		//		/* set to just output and quit */
		//		improvebehave.outputandquit = 1;
		//	}

		//	/* check if this is s */
		//	if (strcmp(argv[i], "-s") == 0)
		//	{
		//		/* parse and set the options from this file */
		//		parseimproveconfig(argv[i+1], b);
		//	}
		//}
	}


	void TetraMeshImprover::parseimprovestateline(ImproveStats *stats)
	{
		/* smoothing stats */
		stats->nonsmoothattempts = 0;
		stats->nonsmoothsuccesses = 0;
		stats->freesmoothattempts = 0;
		stats->freesmoothsuccesses = 0;
		stats->facetsmoothattempts = 0;
		stats->facetsmoothsuccesses = 0;
		stats->segmentsmoothattempts = 0;
		stats->segmentsmoothsuccesses = 0;
		stats->fixedsmoothattempts = 0;
		stats->fixedsmoothsuccesses = 0;

		/* topological stats */
		stats->edgeremovals = 0;
		stats->boundaryedgeremovals = 0;
		stats->edgeremovalattempts = 0;
		stats->boundaryedgeremovalattempts = 0;
		for (int i = 0; i < MAXRINGTETS; i++)
		{
			stats->ringsizesuccess[i] = 0;
			stats->ringsizeattempts[i] = 0;
		}
		stats->faceremovals = 0;
		stats->faceremovalattempts = 0;
		for (int i = 0; i < MAXFACETREESIZE; i++)
		{
			stats->facesizesuccess[i] = 0;
			stats->facesizeattempts[i] = 0;
		}
		stats->flip22attempts = 0;
		stats->flip22successes = 0;
		stats->edgecontractionattempts = 0;
		stats->edgecontractions = 0;
		for (int i = 0; i < NUMEDGECASES+1; i++)
		{
			stats->edgecontractcaseatt[i] = 0;
			stats->edgecontractcasesuc[i] = 0;
		}
		for (int i = 0; i < MAXRINGTETS; i++)
		{
			stats->edgecontractringatt[i] = 0;
			stats->edgecontractringsuc[i] = 0;
		}

		/* vertex insertion stats */
		stats->specialinserttempts = 0;
		stats->specialinsertsuccess = 0;
		stats->inserttempts = 0;
		stats->insertsuccess = 0;

		//		/* timing stats */
		//#ifndef NO_TIMER
		//		struct timeval starttime;
		//#endif /* not NO_TIMER */
		stats->totalmsec = 0;
		stats->smoothmsec = 0;
		stats->topomsec = 0;
		stats->contractmsec = 0;
		stats->insertmsec = 0;
		stats->smoothlocalmsec = 0;
		stats->topolocalmsec = 0;
		stats->insertlocalmsec = 0;
		stats->contractlocalmsec = 0;
		stats->biggestcavityusec = 0;
		stats->finalcavityusec = 0;
		stats->cavityimproveusec = 0;

		/* general stats */
		stats->startnumtets = 0;
		stats->finishnumtets = 0;
		stats->startnumverts = 0;
		stats->finishnumverts = 0;
		stats->dynchangedvol = 0;

		/* quality stats */
		stats->finishworstqual = 0;
		stats->startworstqual = 0;
		stats->startminangle = 0;
		stats->startmaxangle = 0;
		for (int i = 0; i < NUMMEANTHRESHOLDS; i++)
		{
			stats->startmeanquals[i] = 0;
			stats->finishmeanquals[i] = 0;
		}
	}


	/* tetrahedral mesh improvement */
	bool TetraMeshImprover::meshImproving(TetraMesh *mesh_, double &qualafter)
	{
		if (mesh_)
			mesh = mesh_;

		// check mesh
		if (!mesh)
		{
			std::cerr << "TETRAMESH IMPROVEMENT : NO MESH !" << std::endl;
			return false;
		}

		std::vector<TetraHandle> tetstack;      /* stack of tets to be improved */
		int passnum = 1;                        /* current improvement pass */
		int roundsnoimprovement = 0;            /* number of passes since mesh improved */
		bool meansuccess = false;               /* thresholded mean success */
		bool minsuccess = false;                /* minimum quality success */
		double bestmeans[NUMMEANTHRESHOLDS];    /* current best thresholded means */
		bool stop, stop1 = false;               /* whether to continue improvement */
		int numdesperate = 0;

//#ifndef NO_TIMER
//		/* timing vars */
//		struct timeval tv0, tv2;
//		struct timezone tz;
//		/* get initial time */
//		gettimeofday(&tv0, &tz);
//		stats.starttime = tv0;
//#endif /* not NO_TIMER */

		/* perform improvement initialization */
		improveinit(mesh, bestmeans);

		//rebuildmesh(mesh);

		/********** INITIAL SMOOTHING AND TOPO PASSES **********/

		/* initial global special vertex insertion pass */
		//stop = pass(TOPOPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);
		//stop = pass(SPECINSERTPASS, mesh, MININSERTIONIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
		/* clean mesh */
		//mesh->clean_garbage();

		/* initial global smoothing pass */
		//stop = pass(SMOOTHPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);

		/* initial global topological improvement pass */
		//stop = pass(TOPOPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);
		/* clean mesh */
		//mesh->clean_garbage();

		/* initial global contraction improvement pass */
		//if (improvebehave.anisotropic == false)
		//{
		//	//stop = pass(CONTRACTALLPASS, mesh, MINCONTRACTIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
		//	/* clean mesh */
		//	//mesh->clean_garbage();
		//}

		/***************** SIZE CONTROL *************************/
		//if (improvebehave.sizing)
		//{
		//	passnum = sizecontrol(mesh, behave, in, vertexpool, argc, argv);
		//}

		/*************** MAIN IMPROVEMENT LOOP ******************/
		while (roundsnoimprovement < STATICMAXPASSES)
		{
			/* if the mesh is already fine, stop improvement */
			//if (stop) break;

			/* clean mesh */
			//mesh->clean_garbage();
			//rebuildmesh(mesh);

			/* initial global special vertex insertion pass */
			//stop = pass(SPECINSERTPASS, mesh, MININSERTIONIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
			stop = pass(TOPOPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);


			stop = pass(CONTRACTPASS, mesh, MINCONTRACTIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);

			/* perform a smoothing pass */
			stop = pass(SMOOTHPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);

			///* if the smoothing pass failed to sufficiently improve the mesh */
			//if ((minsuccess == false)/* && (meansuccess == false)*/)
			//{
			//	mesh->clean_garbage();

			//	/* perform a global topological pass */
			//	stop = pass(TOPOPASS, mesh, MINSMOOTHITERIMPROVE, minsuccess, meansuccess, passnum++, bestmeans);
			//	/* clean mesh */
			//	//mesh->clean_garbage();
			//	
			//	//stop = pass(SPECINSERTPASS, mesh, MININSERTIONIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);

			//	if (stop) break;

			//	/* if the topo pass also failed */
			//	//if ((minsuccess == false)/* && (meansuccess == false)*/)
			//	//{
			//	//	/* perform a contraction and insertion pass */
			//	//	if (improvebehave.enableinsert)
			//	//	{
			//	//		/* potentially start with a pass of edge contraction */
			//	//		//if (improvebehave.edgecontraction)
			//	//		//{
			//	//		//	stop = pass(CONTRACTPASS, mesh, MINCONTRACTIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
			//	//		//	/* clean mesh */
			//	//		//	//mesh->clean_garbage();
			//	//		//}
			//	//		//else
			//	//		//{
			//	//		//	stop = false;
			//	//		//}

			//	//		//if (roundsnoimprovement == 1 && numdesperate < DESPERATEMAXPASSES)
			//	//		//{
			//	//		//	stop = pass(DESPERATEPASS, mesh, MINCONTRACTIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
			//	//		//	numdesperate++;
			//	//		//	/* clean mesh */
			//	//		//	//mesh->clean_garbage();
			//	//		//	if (stop || stop1) break;
			//	//		//}
			//	//		//else
			//	//		//{
			//	//			//stop = pass(INSERTPASS, mesh, MININSERTIONIMPROVEMENT, minsuccess, meansuccess, passnum++, bestmeans);
			//	//			//if (stop || stop1) break;
			//	//		//}
			//	//	}
			//	//}
			//}
			/* if this pass failed to see any improvement, note it */
			roundsnoimprovement++;
			//if ((minsuccess == false) /*&& (meansuccess == false)*/)
			//{
			//	roundsnoimprovement++;
			//}
			///* reset number of rounds on smoothing success */
			//else 
			//{
			//	roundsnoimprovement = 0;
			//}
		}

		/******************** END MAIN IMPROVEMENT LOOP **********************/

//#ifndef NO_TIMER
//		/* record total time */
//		gettimeofday(&tv2, &tz);
//		stats.totalmsec = msecelapsed(tv0, tv2);
//#endif /* not NO_TIMER */
//
//		/* perform post-improvement cleanup */
		//improvedeinit(mesh, vertexpool, &tetstack, behave, in, argc, argv);
		//mesh->clean_garbage();

		qualafter = qualcalculator.minmeshquality(mesh, qualitymetric);

		return true;
	}


	/* pre-improvement initialization code */
	void TetraMeshImprover::improveinit(TetraMesh *mesh, double bestmeans[NUMMEANTHRESHOLDS])
	{
		//int consistent;
		//double minqualbefore;
		double meanqualbefore[NUMMEANTHRESHOLDS];
		double worstin;
		double worstinitqual;
		TetraMesh::HedronIter h_iter;
		TetraMesh::PointIter p_iter;

		for (int i = 0; i < NUMMEANTHRESHOLDS; i++)
		{
			bestmeans[i] = 0.0;
		}

		stats.startnumverts = mesh->size_point();
		stats.startnumtets = mesh->size_tetrahedron();

		///* stack of tets to be improved */
		//improveTetraVec.clear();
		//for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		//{
		//	improveTetraVec.push_back(h_iter.handle());
		//}

		///* stack of points to be improved */
		//improvePointVec.clear();
		//for (p_iter = mesh->points_begin(); p_iter != mesh->points_end(); ++p_iter)
		//{
		//	improvePointVec.push_back(p_iter.handle());
		//}
		

		/* this stack stores a journal of improvement steps */
		journals.clear();
		journalentry = 0;

		/* get the worst input angle */
		worstin = qualcalculator.worstinputangle(mesh);

		/* get the worst input quality */
		worstinitqual = qualcalculator.worstquality(mesh);

		/* build stack for initial quality evaluation */
		meanimprove(bestmeans, meanqualbefore, SMOOTHPASS);

		/* set initial minimum and thresholded mean qualities */
		for (int i = 0; i < NUMMEANTHRESHOLDS; i++)
		{
			stats.startmeanquals[i] = meanqualbefore[i];
		}
		stats.startworstqual = worstinitqual;
	}

	/* see if new means contains any better means. if so,
	update best means and return true */
	bool TetraMeshImprover::meanimprove(double bestmeans[], double newmeans[], int passtype)
	{
		int i;
		bool foundbetter = false;
		double minimprovement = improvebehave.minstepimprovement;

		if (passtype == INSERTPASS)
		{
			minimprovement = improvebehave.mininsertionimprovement;
		}
		if (passtype == DESPERATEPASS)
		{
			minimprovement = improvebehave.mininsertionimprovement;
		}

		for (i=0; i<NUMMEANTHRESHOLDS; i++)
		{
			if (newmeans[i] > bestmeans[i])
			{
				/* see if it beats it by the required threshold */
				if (newmeans[i] - bestmeans[i] > minimprovement)
					foundbetter = true;
				bestmeans[i] = newmeans[i];
			}
		}
		return foundbetter;
	}

	/* run a pass (smoothing, topo, insertion). return true
	if we have reached the desired quality */
	bool TetraMeshImprover::pass(int passtype, TetraMesh *mesh, double threshold, 
		                         bool &minsuccess, bool &meansuccess,
		                         int passnum, double bestmeans[])
	{
		/* quality vars */
		double minqualbefore, minqualafter;
		double meanqualbefore[NUMMEANTHRESHOLDS], meanqualafter[NUMMEANTHRESHOLDS];
		//double minedge, maxedge, meanedge;

		/* smoothing vars */
		int smoothkinds = 0;

		/* sine of target angles */
		double goalangleminsine = degtosin(improvebehave.goalanglemin);
		double goalanglemaxsine = degtosin(improvebehave.goalanglemax);
		double biggestangle;
		double smallestangle;

		bool desperate = false;
		int passstartid = 0;

		if (improvebehave.facetsmooth) smoothkinds |= SMOOTHFACETVERTICES;
		if (improvebehave.segmentsmooth) smoothkinds |= SMOOTHSEGMENTVERTICES;
		if (improvebehave.fixedsmooth) smoothkinds |= SMOOTHFIXEDVERTICES;

		assert(passtype == SMOOTHPASS ||
			   passtype == TOPOPASS ||
			   passtype == CONTRACTPASS ||
			   passtype == CONTRACTALLPASS ||
			   passtype == INSERTPASS ||
			   passtype == DESPERATEPASS ||
			   passtype == SPECINSERTPASS);

		/* get global worst quality; this is what must improve */
		minqualbefore = qualcalculator.minmeshquality(mesh, qualitymetric);

		/* capture animation info */
		if (improvebehave.animate)
		{
			passstartid = journals.size(); 
		}

		/* run the actual pass */
		switch (passtype)
		{
		case SMOOTHPASS:
			//vnonsmoother.initialinformation();
			//vnonsmoother.setQualityKind(MINSINE2);
			//vnonsmoother.setBoundaryChange(false);
			//vnonsmoother.setSmoothingType(0x11);
			//vnonsmoother.NonsmoothVertex(1.0e-6);
			//vnonsmoother.upDate();
			//minqualafter = qualcalculator.minmeshquality(mesh, qualitymetric);
			vsmoother.smoothpass(NULL, threshold, bestmeans, 
				meanqualafter, &minqualafter, smoothkinds);
			break;
		case TOPOPASS:
			toperator.topologyPass(NULL, NULL, minqualafter);
			break;
		case CONTRACTPASS:
			toperator.contractworst(improvebehave.insertthreshold, bestmeans, 
				                    meanqualafter, &minqualafter, true);
			break;
		case CONTRACTALLPASS:
			toperator.contractpass(NULL, NULL, minqualafter, true);
			break;
		case DESPERATEPASS:
		case INSERTPASS:
			vinserter.targetTetras(minqualafter);
			break;
		default:
			printf("i don't know how to run pass type %d, dying\n", passtype);
			exit(1);
		}

		/* check for success */
		//meansuccess = meanimprove(bestmeans, meanqualafter, passtype);
		if (minqualafter - minqualbefore < MINMINIMPROVEMENT)
		{
			minsuccess = false;
		}
		else
		{
			minsuccess = true;
		}

		/* check whether we have reached the goal quality for minimum or maximum angle */
		if (minqualafter > goalangleminsine || minqualafter > goalanglemaxsine)
		{
			/* compute the extreme angles */
			getextremeangles(mesh, smallestangle, biggestangle);

			/* we must have reached one of these angles */
			/* not necessarily true for non-angle based quality measures */
			/* assert(smallestangle > improvebehave.goalanglemin || biggestangle < improvebehave.goalanglemax); */

			/* if we've reached both thresholds, let the main loop know we're done */
			if (smallestangle > improvebehave.goalanglemin && biggestangle < improvebehave.goalanglemax)
			{
				return true;
			}
		}
		return false;
	}


	void TetraMeshImprover::getextremeangles(TetraMesh *mesh, double &outsmallestangle, double &outbiggestangle, 
		                                     int stat[180])
	{
		Point points[4];
		TetraMesh::Normal t, u, v;
		TetraMesh::Normal facenormals[4];
		double tetangles[6];
		double smallangle, bigangle;
		TetraMesh::HedronIter h_it;

		for (int i = 0; i < 180; i++)
		{
			stat[i] = 0;
		}
		// for each tetras calculate its smallest and biggest angles
		smallangle = 180.0;
		bigangle = 0.0;
		for (h_it = mesh->hedrons_begin(); h_it != mesh->hedrons_end(); ++h_it)
		{
			// skip invalid tetrahedrons
			if (!mesh->is_valid(h_it.handle()))
				continue;

			// fetch tetra points
			tetPoints(mesh, h_it.handle(), points);
			// get t, u, v vectors
			t = points[1] - points[0];
			u = points[2] - points[0];
			v = points[3] - points[0];
			// get four normalized face normals
			facenormals[0] = ((u - t) % (v - t)).normalize();
			facenormals[1] = (v % u).normalize();
			facenormals[2] = (t % v).normalize();
			facenormals[3] = (u % t).normalize();
			// get six tetra dihedral angles
			tetangles[0] = 180.0 - acos(facenormals[0] | facenormals[1]) * 180.0 / PI;
			tetangles[1] = 180.0 - acos(facenormals[0] | facenormals[2]) * 180.0 / PI;
			tetangles[2] = 180.0 - acos(facenormals[0] | facenormals[3]) * 180.0 / PI;
			tetangles[3] = 180.0 - acos(facenormals[1] | facenormals[2]) * 180.0 / PI;
			tetangles[4] = 180.0 - acos(facenormals[1] | facenormals[3]) * 180.0 / PI;
			tetangles[5] = 180.0 - acos(facenormals[2] | facenormals[3]) * 180.0 / PI;

			// fetch statistics of angle and smallest angle and biggest angle
			for (int i = 0; i < 6; i ++)
			{
				if (tetangles[i] < smallangle)
					smallangle = tetangles[i];
				if (tetangles[i] > bigangle)
					bigangle = tetangles[i];

				for (int j = 1; j < 181; j++)
				{
					if (tetangles[i] < j)
					{
						++ stat[j-1];
						break;
					}
				}
			}
		}
		outsmallestangle = smallangle;
		outbiggestangle = bigangle;
	}

	void TetraMeshImprover::getextremeangles(TetraMesh *mesh, double &outsmallestangle, double &outbiggestangle)
	{
		Point points[4];
		TetraMesh::Normal t, u, v;
		TetraMesh::Normal facenormals[4];
		double tetangles[6];
		double smallangle, bigangle;
		TetraMesh::HedronIter h_it;

		// for each tetras calculate its smallest and biggest angles
		smallangle = 180.0;
		bigangle = 0.0;
		for (h_it = mesh->hedrons_begin(); h_it != mesh->hedrons_end(); ++h_it)
		{
			// skip invalid tetrahedrons
			if (!mesh->is_valid(h_it.handle()))
				continue;

			// fetch tetra points
			tetPoints(mesh, h_it.handle(), points);
			// get t, u, v vectors
			t = points[1] - points[0];
			u = points[2] - points[0];
			v = points[3] - points[0];
			// get four normalized face normals
			facenormals[0] = ((u - t) % (v - t)).normalize();
			facenormals[1] = (v % u).normalize();
			facenormals[2] = (t % v).normalize();
			facenormals[3] = (u % t).normalize();
			// get six tetra dihedral angles
			tetangles[0] = 180.0 - acos(facenormals[0] | facenormals[1]) * 180.0 / PI;
			tetangles[1] = 180.0 - acos(facenormals[0] | facenormals[2]) * 180.0 / PI;
			tetangles[2] = 180.0 - acos(facenormals[0] | facenormals[3]) * 180.0 / PI;
			tetangles[3] = 180.0 - acos(facenormals[1] | facenormals[2]) * 180.0 / PI;
			tetangles[4] = 180.0 - acos(facenormals[1] | facenormals[3]) * 180.0 / PI;
			tetangles[5] = 180.0 - acos(facenormals[2] | facenormals[3]) * 180.0 / PI;

			// fetch statistics of angle and smallest angle and biggest angle
			for (int i = 0; i < 6; i ++)
			{
				if (tetangles[i] < smallangle)
					smallangle = tetangles[i];
				if (tetangles[i] > bigangle)
					bigangle = tetangles[i];
			}
		}

		outsmallestangle = smallangle;
		outbiggestangle = bigangle;
	}

	//void TetraMeshImprover::getextremeangles(TetraMesh *mesh, double &outsmallestangle, double &outbiggestangle)
	//{
	//	Point point[4];
	//	Point tcenter;
	//	double tansquaretable[8];
	//	double aspecttable[16];
	//	double circumtable[16];
	//	double edgelength[3][4];
	//	double dx, dy, dz;
	//	double dotproduct;
	//	double tansquare;
	//	double pyrvolume;
	//	double facearea2;
	//	double pyrbiggestface2;
	//	double shortest, longest;
	//	double pyrshortest2, pyrlongest2;
	//	double smallestvolume, biggestvolume;
	//	double pyrminaltitude2;
	//	double minaltitude;
	//	double pyraspect2;
	//	double worstaspect;
	//	double pyrcircumratio2;
	//	double worstcircumratio;
	//	double smallestangle, biggestangle;
	//	double radconst, degconst;
	//	int anglecount[18];
	//	int aspectcount[16];
	//	int circumcount[16];
	//	int aspectindex;
	//	int circumindex;
	//	int tendegree;
	//	int acutebiggestflag;
	//	int firstiterationflag;
	//	int i, ii, j, k, l;
	//	TetraMesh::Normal E[3];
	//	TetraMesh::Normal facenormal[4];
	//	TetraMesh::HedronIter h_iter;

	//	radconst = PI / 18.0;
	//	degconst = 180.0 / PI;
	//	for (i = 0; i < 8; i++) 
	//	{
	//		tansquaretable[i] = tan(radconst * (double) (i + 1));
	//		tansquaretable[i] = tansquaretable[i] * tansquaretable[i];
	//	}
	//	for (i = 0; i < 18; i++) 
	//	{
	//		anglecount[i] = 0;
	//	}

	//	aspecttable[0]  =      1.5;      aspecttable[1]  =     2.0;
	//	aspecttable[2]  =      2.5;      aspecttable[3]  =     3.0;
	//	aspecttable[4]  =      4.0;      aspecttable[5]  =     6.0;
	//	aspecttable[6]  =     10.0;      aspecttable[7]  =    15.0;
	//	aspecttable[8]  =     25.0;      aspecttable[9]  =    50.0;
	//	aspecttable[10] =    100.0;      aspecttable[11] =   300.0;
	//	aspecttable[12] =   1000.0;      aspecttable[13] = 10000.0;
	//	aspecttable[14] = 100000.0;      aspecttable[15] =     0.0;
	//	for (i = 0; i < 16; i++) {
	//		aspectcount[i] = 0;
	//	}

	//	circumtable[0]  =      0.75;     circumtable[1]  =     1.0;
	//	circumtable[2]  =      1.5;      circumtable[3]  =     2.0;
	//	circumtable[4]  =      3.0;      circumtable[5]  =     5.0;
	//	circumtable[6]  =     10.0;      circumtable[7]  =    15.0;
	//	circumtable[8]  =     25.0;      circumtable[9]  =    50.0;
	//	circumtable[10] =    100.0;      circumtable[11] =   300.0;
	//	circumtable[12] =   1000.0;      circumtable[13] = 10000.0;
	//	circumtable[14] = 100000.0;      circumtable[15] =     0.0;
	//	for (i = 0; i < 16; i++) {
	//		circumcount[i] = 0;
	//	}

	//	shortest = 0.0;
	//	longest = 0.0;
	//	smallestvolume = 0.0;
	//	biggestvolume = 0.0;
	//	minaltitude = 0.0;
	//	worstaspect = 0.0;
	//	worstcircumratio = 0.0;
	//	smallestangle = 100.0;
	//	biggestangle = 0.0;
	//	acutebiggestflag = 1;

	//	firstiterationflag = 1;
	//	for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
	//	{
	//		if (!mesh->is_valid(h_iter.handle()))
	//			continue;

	//		tetPoints(mesh, h_iter.handle(), point);

	//		pyrshortest2 = 0.0;
	//		pyrlongest2 = 0.0;
	//		pyrbiggestface2 = 0.0;

	//		for (i = 0; i < 4; i++) {
	//			j = (i + 1) & 3;
	//			if ((i & 1) == 0) {
	//				k = (i + 3) & 3;
	//				l = (i + 2) & 3;
	//			} else {
	//				k = (i + 2) & 3;
	//				l = (i + 3) & 3;
	//			}

	//			facenormal[i][0] =
	//				(point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
	//				(point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
	//			facenormal[i][1] =
	//				(point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
	//				(point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
	//			facenormal[i][2] =
	//				(point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
	//				(point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
	//			facearea2 = facenormal[i][0] * facenormal[i][0] +
	//				facenormal[i][1] * facenormal[i][1] +
	//				facenormal[i][2] * facenormal[i][2];
	//			if (facearea2 > pyrbiggestface2) {
	//				pyrbiggestface2 = facearea2;
	//			}

	//			for (j = i + 1; j < 4; j++) {
	//				dx = point[i][0] - point[j][0];
	//				dy = point[i][1] - point[j][1];
	//				dz = point[i][2] - point[j][2];
	//				edgelength[i][j] = dx * dx + dy * dy + dz * dz;
	//				if (edgelength[i][j] > pyrlongest2) {
	//					pyrlongest2 = edgelength[i][j];
	//				}
	//				if ((j == 1) || (edgelength[i][j] < pyrshortest2)) {
	//					pyrshortest2 = edgelength[i][j];
	//				}
	//			}
	//		}
	//		if (pyrlongest2 > longest) {
	//			longest = pyrlongest2;
	//		}
	//		if ((pyrshortest2 < shortest) || firstiterationflag) {
	//			shortest = pyrshortest2;
	//		}

	//		pyrvolume = (double)orient(point[0], point[1], point[2], point[3]);
	//		if ((pyrvolume < smallestvolume) || firstiterationflag) {
	//			smallestvolume = pyrvolume;
	//		}
	//		if (pyrvolume > biggestvolume) {
	//			biggestvolume = pyrvolume;
	//		}
	//		pyrminaltitude2 = pyrvolume * pyrvolume / pyrbiggestface2;
	//		if ((pyrminaltitude2 < minaltitude) || firstiterationflag) {
	//			minaltitude = pyrminaltitude2;
	//		}
	//		pyraspect2 = pyrlongest2 / pyrminaltitude2;
	//		if (pyraspect2 > worstaspect) {
	//			worstaspect = pyraspect2;
	//		}
	//		aspectindex = 0;
	//		while ((pyraspect2 >
	//			aspecttable[aspectindex] * aspecttable[aspectindex]) &&
	//			(aspectindex < 15)) {
	//				aspectindex++;
	//		}
	//		aspectcount[aspectindex]++;

	//		tcenter = (point[0] + point[1] + point[2] + point[3])/4.0;
	//		pyrcircumratio2 = (tcenter[0] * tcenter[0] + tcenter[1] * tcenter[1] +
	//			tcenter[2] * tcenter[2]) / pyrshortest2;
	//		if (pyrcircumratio2 > worstcircumratio) {
	//			worstcircumratio = pyrcircumratio2;
	//		}
	//		circumindex = 0;
	//		while ((pyrcircumratio2 >
	//			circumtable[circumindex] * circumtable[circumindex]) &&
	//			(circumindex < 15)) {
	//				circumindex++;
	//		}
	//		circumcount[circumindex]++;

	//		for (i = 0; i < 3; i++) {
	//			for (j = i + 1; j < 4; j++) {
	//				k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
	//				l = 6 - i - j - k;
	//				dotproduct = facenormal[i][0] * facenormal[j][0] +
	//					facenormal[i][1] * facenormal[j][1] +
	//					facenormal[i][2] * facenormal[j][2];
	//				if (dotproduct == 0.0) {
	//					anglecount[9]++;
	//					if (acutebiggestflag)
	//					{
	//						biggestangle = 1.0e+300;
	//						acutebiggestflag = 0;
	//					}
	//				} else {
	//					tansquare = pyrvolume * pyrvolume * edgelength[k][l] /
	//						(dotproduct * dotproduct);
	//					tendegree = 8;
	//					for (ii = 7; ii >= 0; ii--) {
	//						if (tansquare < tansquaretable[ii]) {
	//							tendegree = ii;
	//						}
	//					}
	//					if (dotproduct < 0.0) {
	//						anglecount[tendegree]++;
	//						if (tansquare < smallestangle) {
	//							smallestangle = tansquare;
	//						}
	//						if (acutebiggestflag && (tansquare > biggestangle)) {
	//							biggestangle = tansquare;
	//						}
	//					} else {
	//						anglecount[17 - tendegree]++;
	//						if (acutebiggestflag || (tansquare < biggestangle)) {
	//							biggestangle = tansquare;
	//							acutebiggestflag = 0;
	//						}
	//					}
	//				}
	//			}
	//		}
	//		firstiterationflag = 0;
	//	}

	//	shortest = sqrt(shortest);
	//	longest = sqrt(longest);
	//	minaltitude = sqrt(minaltitude);
	//	worstaspect = sqrt(worstaspect);
	//	worstcircumratio = sqrt(worstcircumratio);
	//	smallestvolume /= 6.0;
	//	biggestvolume /= 6.0;
	//	smallestangle = degconst * atan(sqrt(smallestangle));
	//	if (acutebiggestflag) {
	//		biggestangle = degconst * atan(sqrt(biggestangle));
	//	} else {
	//		biggestangle = 180.0 - degconst * atan(sqrt(biggestangle));
	//	}
	//	outsmallestangle = smallestangle;
	//	outbiggestangle = biggestangle;
	//}


	// rebuild mesh
	void TetraMeshImprover::rebuildmesh(TetraMesh *mesh)
	{
		TetraMesh *oldmesh;            // tetramesh pointer point to the old mesh
		PContainer pc;
		VContainer vc;
		TetraMesh::PointIter p_it;
		TetraMesh::HedronIter h_it;
		TetraMesh::HedronVertexIter hv_it;

		// save old mesh
		oldmesh = mesh;

		// clean old mesh
		oldmesh->clean_garbage();

		// create new mesh
		mesh = new TetraMesh;
		this->mesh = mesh;

		// fetch mesh points
		for (p_it = oldmesh->points_begin(); p_it != oldmesh->points_end(); ++p_it)
		{
			pc.push_back(oldmesh->point(p_it.handle()));
		}

		// fetch hedron vertex
		for (h_it = oldmesh->hedrons_begin(); h_it != oldmesh->hedrons_end(); ++h_it)
		{
			for (hv_it = oldmesh->hedron_vertex_iter(h_it.handle()); hv_it; ++hv_it)
			{
				vc.push_back(Vertex(hv_it.handle(),oldmesh->point_handle(hv_it.handle())));
			}
		}

		// build new mesh
		mesh->set_point_container(pc);
		mesh->set_vertex_container(vc);
		mesh->build_topology();
		mesh->request_face_normals();
		mesh->update_face_normals();

		// set tools' mesh
		vsmoother.setMesh(mesh);
		vinserter.setMesh(mesh);
		quadric.setMesh(mesh);
		toperator.setMesh(mesh);
		qualcalculator.setMesh(mesh);

		// delete old mesh
		delete oldmesh;

	}
}