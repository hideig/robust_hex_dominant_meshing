#ifndef _TOPOLOGICAL_OPERATOR_TOP_DATA_DEFINE
#define _TOPOLOGICAL_OPERATOR_TOP_DATA_DEFINE

#include <VolumeMesh/Mesh/TetraMesh.h>
#include <iterator>
#include <deque>
#include <algorithm>
namespace VolumeMesh
{

	/*****************************************************************************/
	/*                                                                           */
	/*  Defines                                                                  */
	/*                                                                           */
	/*****************************************************************************/

	/* a really big floating point number */
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

	/* do paranoid checks during improvement */
#define IMPROVEPARANOID false

	/* verbosity of debugging output */
#define IMPROVEVERBOSITY 5

	/* a really big floating point number */
#define HUGEFLOAT 1.0e100

	/* the sines of some angles */
#define SINEHALF 0.00872653549837
#define SINE1  0.01745240644
#define SINE2  0.0348994967
#define SINE5  0.08715574275
#define SINE10 0.17364817767
#define SINE13 0.22495105434
#define SINE15 0.2588190451
#define SINE20 0.34202014333
#define SINE25 0.42261826174
#define SINE30 0.5
#define SINE35 0.57357643635
#define SINE40 0.64278760969
#define SINE45 0.70710678119
#define SINEEQUILATERAL 0.94280903946

	/* when the journal reaches this size, half it's size (remove the older half
	of the entries) */
#define JOURNALHALFSIZE 1000000

	/* vertex types */
#define INPUTVERTEX 0
#define FIXEDVERTEX 1
#define SEGMENTVERTEX 2
#define FACETVERTEX 3
#define FREEVERTEX 4
#define UNDEADVERTEX 15
	/* new kind of vertex to identify the ones I've put in */
#define INSERTEDVERTEX 4

	/* types of improvement passes */
#define SMOOTHPASS 0
#define TOPOPASS 1
#define CONTRACTPASS 2
#define INSERTPASS 3
#define DESPERATEPASS 4
#define DEFORMPASS 5
#define SPECINSERTPASS 6
	/* types of local improvement passes */
#define SMOOTHPASSLOCAL 6
#define TOPOPASSLOCAL 7
#define CONTRACTPASSLOCAL 8
#define INSERTPASSLOCAL 9
	/* size control pass */
#define SIZECONTROLPASS 10
#define CONTRACTALLPASS 11

	/* number of quality measures */
#define NUMQUALMEASURES 6

	/* number of "thresholded" means used to approximate mesh quality */
#define NUMMEANTHRESHOLDS 7

	/* edge cases */
#define NUMEDGECASES 10
#define NOEDGECASE 0
#define FREEFREEEDGE 1
#define FREEFACETEDGE 2
#define FREESEGMENTEDGE 3
#define FREEFIXEDEDGE 4
#define FACETFACETEDGE 5
#define FACETSEGMENTEDGE 6
#define FACETFIXEDEDGE 7
#define SEGMENTSEGMENTEDGE 8
#define SEGMENTFIXEDEDGE 9
#define FIXEDFIXEDEDGE 10

	/*vertex insert*/
#define EDGELABEL 0
//#define TETLABEL 1
#define CAVLABEL 0
#define ANTICAVLABEL 1
#define NOLABEL 2
#define DEPTHTABLESIZE 10
#define NOCAVITYTET TetraHandle(-1)
#define NOCAVITYFACE TetraMesh::HalfFaceHandle(-1)
#define GHOSTTET TetraHandle(-1)
#define EPSILON 1E-5
#define O3DERRORBOUNDA (7.0 + 56.0 * EPSILON) * EPSILON


	/* number of passes without improvement before static improvement quits */
#define STATICMAXPASSES 5
	/* number of desperate insertion passes that can ever be attempted */
#define DESPERATEMAXPASSES 3

	/*****************************************************************************/
	/*  Topological improvement options                                          */
	/*****************************************************************************/

	/* number of tets to allow in a ring around an edge to be removed */
#define MAXRINGTETS 70
#define MAXRINGTETS2 50

	/* maximum number of tets in sandwich set replacement during edge removal */
#define MAXNEWTETS 150

	/* maximum number of faces in tree for multi-face removal */
#define MAXFACETREESIZE 50

	/* minimum quality to allow a 4-verts-on-boundary tet to be created
	by a topological improvement operation */
#define MIN4BOUNDQUAL SINE1

	/*****************************************************************************/
	/*  Smoothing improvement options                                            */
	/*****************************************************************************/

	/* do paranoid checks of smoothing process */
#define SMOOTHPARANOID false

	/* how close to the worst does a function have to be to be in the
	active set? */
#define ACTIVESETFACTOR 1.03

	/* how many tets will we allow incident to one vertex for smoothing? */
#define MAXINCIDENTTETS 700

	/* minimum step size */
#define MINSTEPSIZE 1.0e-5

	/* maximum iterations in non-smooth line search */
#define MAXLINEITER 50

	/* rate must be worse by this much to reset alpha */
#define RATEEPSILON 1.0e-6

	/* maximum iterations of non smooth optimization */
#define MAXSMOOTHITER 50

	/* minimum quality improvement in a smoothing step */
#define MINSMOOTHITERIMPROVE 1.0e-5

	/* if d is within this of zero-length, call it zero */
#define DEPSILON 1.0e-5

	/* if closest point computation has a factor smaller than this, make it the origin */
#define NEARESTMIN 1.0e-13

	/* the minimum improvement of the minimum quality element for a 
	(topo|smoothing|insertion) pass to "succeed" */
#define MINMINIMPROVEMENT 1.0e-6

#define MINLOCALMEANIMPROVE 0.005

	/* minimum improvement required for edge contraction to succeed */
#define MINCONTRACTIMPROVEMENT 1.0e-06

	/* the dot product of two normals must be 1+- this value to be considered coplanar */
#define COPLANARTOL 1e-4

	/* the absolute value of the dot product of two edges must be 1+- 
	this value to be considered colinear */
#define COLINEARTOL 1e-4

	/* determines whether facet/segment vertices are smoothed */
#define SMOOTHFACETVERTICES 0x01
#define SMOOTHSEGMENTVERTICES 0x02
#define SMOOTHFIXEDVERTICES 0x04

	/*****************************************************************************/
	/*  Insertion/cavity improvement options                                     */
	/*****************************************************************************/

	/* do paranoid checks during insertion */
#define INSERTPARANOID false

	/* minimum improvement for a submesh */
#define MINSUBMESHIMPROVEMENT 1.0e-3

	/* maximum submesh improvement iterations */
#define MAXSUBMESHITERATIONS 8

	/* if a tet is within this number of the worst tet, its close to worst */
#define CLOSETOWORST SINE2

	/* how much bigger to allow the improvement stack to be when we've got 
	one of the worst tets */ 
#define TRYHARDFACTOR 15
#define TRYHARDMAXSUBMESHITERATIONS 20
#define TRYHARDMINSUBMESHIMPROVEMENT 1e-10

	/* minimum quality of an intermediate tet */
#define MINTETQUALITY 1.0e-14
#define MINSIZETETQUALITY 1.0e-10

	/* maximum size of stuff for cavity drilling */
#define MAXCAVITYFACES 10000
#define MAXCAVITYTETS 10000

	/* minimum improvement of vertex insertion */
#define MININSERTIONIMPROVEMENT 1.0e-13

	/* minimum positivity of orientation for a face to be oriented "toward" */
#define MINFACING 1.0e-7

	/* factor by which child cavity qualities are increased */
#define CHILDFAVORFACTOR 1.0

	/* maximum difference between two qualities where they are considered 'equal' */
#define MAXQUALDIFFERENCE 1.0e-15

	/* maximum number of outgoing faces allowed for a tet */
#define MAXOUTFACES 200

	/* maximum number of edges in a cavity dag */
#define MAXCAVITYDAGEDGES 50000

	/* perform initial/final smooth of all cavity vertices */
#define INITIALCAVITYSMOOTH true
#define FINALCAVITYSMOOTH true

	/* maximum stack size for cavity improvement */
#define MAXCAVITYSTACK 250

	/* biggest size for cavity heap */ 
#define MAXCAVITYHEAPSIZE MAXCAVITYTETS

	/* deepest level a tet can be */
#define MAXCAVDEPTH 1000

	/* failure count on which to perform desperate insertion pass */
#define DESPERATEPASSNUM 2
	/* if the quality of the mesh is really terrible, don't bother trying to insert 
	up to max insert quality */
#define QUALFROMDESPERATE SINE15
	/* if the minimum quality is high, reach even higher by this much in desperate pass */
#define QUALUPFROMDESPERATE SINE1

	/* never attempt insertion on more than this many tets, no matter what */
#define MAXINSERTTETS 4000

	/*****************************************************************************/
	/*  size control options                                                     */
	/*****************************************************************************/

	/* minimum quality of affect tets after size-control edge contraction */
#define MINCONTRACTQUAL 5.0e-2
#define MINSIZESPLITQUAL 5.0e-2
#define MEANEDGESCALE 0.8
	/* maximum number of size control iterations */
	/*#define MAXSIZEITERS 30*/
#define MAXSIZEITERS 10
#define CONTROLSHORTFAC 1.35
#define CONTROLLONGFAC 0.85

	/*****************************************************************************/
	/*                                                                           */
	/*  Data structures                                                          */
	/*                                                                           */
	/*****************************************************************************/

	/* structure to hold global improvement statistics */
	struct ImproveStats
	{
		/* smoothing stats */
		int nonsmoothattempts;
		int nonsmoothsuccesses;
		int freesmoothattempts;
		int freesmoothsuccesses;
		int facetsmoothattempts;
		int facetsmoothsuccesses;
		int segmentsmoothattempts;
		int segmentsmoothsuccesses;
		int fixedsmoothattempts;
		int fixedsmoothsuccesses;

		/* topological stats */
		int edgeremovals;
		int boundaryedgeremovals;
		int edgeremovalattempts;
		int boundaryedgeremovalattempts;
		int ringsizesuccess[MAXRINGTETS];
		int ringsizeattempts[MAXRINGTETS];
		int faceremovals;
		int faceremovalattempts;
		int facesizesuccess[MAXFACETREESIZE];
		int facesizeattempts[MAXFACETREESIZE];
		int flip22attempts;
		int flip22successes;
		int edgecontractionattempts;
		int edgecontractions;
		int edgecontractcaseatt[NUMEDGECASES+1];
		int edgecontractcasesuc[NUMEDGECASES+1];
		int edgecontractringatt[MAXRINGTETS];
		int edgecontractringsuc[MAXRINGTETS];
		int specialinserttempts;
		int specialinsertsuccess;
		int inserttempts;
		int insertsuccess;

//		/* timing stats */
//#ifndef NO_TIMER
//		struct timeval starttime;
//#endif /* not NO_TIMER */
		int totalmsec;
		int smoothmsec;
		int topomsec;
		int contractmsec;
		int insertmsec;
		int smoothlocalmsec;
		int topolocalmsec;
		int insertlocalmsec;
		int contractlocalmsec;
		int biggestcavityusec;
		int finalcavityusec;
		int cavityimproveusec;

		/* general stats */
		int startnumtets;
		int finishnumtets;
		int startnumverts;
		int finishnumverts;
		double dynchangedvol;

		/* quality stats */
		double finishworstqual;
		double startworstqual;
		double startminangle;
		double startmaxangle;
		double startmeanquals[NUMMEANTHRESHOLDS];
		double finishmeanquals[NUMMEANTHRESHOLDS];
	};

	struct CavityFace
	{
		TetraMesh::HalfFaceHandle handle;
		double qual;
		TetraHandle child;
		bool inH;
	};

	struct CavityTet
	{
		TetraHandle handle;
		double qual;
		int depth;
		std::vector<CavityFace> outfaces;
		TetraHandle parents[3];
		int label;
	};

	struct CavityEdgeorTet
	{
		double qual;
		int label;
		TetraHandle parent;
		TetraHandle child;
		int childnum;
	};

	///* structure to hold all the options for improvement */
	struct ImproveBehavior
	{
		/* Quality measure */
		int qualmeasure;             /* quality measure used */
		double sinewarpfactor;     /* warp obtuse sines by this factor */

		/* Quadric smoothing options */
		int usequadrics;             /* incorporate quadric error into the objective function */
		double quadricoffset;      /* quality to start every quadric at */
		double quadricscale;       /* factor to scale quadric by */

		/* Smoothing options */
		int nonsmooth;               /* enable non-smooth optimization-based smoothing */
		int facetsmooth;             /* enable smoothing of facet vertices */
		int segmentsmooth;           /* enable smoothing of segment vertices */
		int fixedsmooth;             /* enable smoothing of fixed vertices */

		/* Topological options */
		int edgeremoval;             /* enable edge removal */
		int edgecontraction;         /* enable edge contraction */
		int boundedgeremoval;        /* enable boundary edge removal */
		int singlefaceremoval;       /* enable single face removal (2-3 flips) */
		int multifaceremoval;        /* enable multi face removal */
		int flip22;                  /* enable 2-2 flips */
		int jflips;                  /* use Jonathan's faster flip routines */

		/* Insertion options */
		int enableinsert;            /* global enable of insertion */
		double insertthreshold;    /* percent worst tets */
		int insertbody;              /* enable body vertex insertion */
		int insertfacet;             /* enable facet insertion */
		int insertsegment;           /* enablem segment insertion */
		int cavityconsiderdeleted;   /* consider enlarging cavity for deleted tets? */
		int cavdepthlimit;           /* only allow initial cavity to includes tets this deep */

		/* Special Insertion options */
		int enablespecialinsert;
		double specialinsertthreshold;
		int specialinsertbody;
		int specialinsertfacet;

		/* anisotropic meshing options */
		int anisotropic;             /* globally enable space warping with deformation tensor */
		int tensor;                  /* which scaling tensor field to use */
		int tensorb;                 /* second tensor, to blend with the first one */
		double tensorblend;  /* between 0 and 1, how much of anisotropy to use */

		/* sizing options */
		int sizing;                  /* globally enable mesh element size control */
		int sizingpass;              /* enable or disable initial edge length control */
		double targetedgelength;   /* edge length of the ideal edge for this mesh */
		double longerfactor;       /* factor by which an edge can be longer */
		double shorterfactor;      /* factor by which an edge can be shorter */

		/* Dynamic improvement options */
		double dynminqual;         /* minimum quality demanded for dynamic improvement. */
		double dynimproveto;       /* after minimum quality is reached, improve to at least this level */
		int deformtype;              /* which fake deformation to use */
		int dynimprove;              /* perform dynamic improvement with fake deformation? */

		/* thresholds */
		double minstepimprovement; /* demand at least this much improvement in the mean per step */
		double mininsertionimprovement; /* demand in improvement for insertion */
		double maxinsertquality[NUMQUALMEASURES];   /* never attempt insertion in a tet better than this */

		/* improvement limits */
		double goalanglemin;       /* stop improvement if smallest angle reaches this threshold */
		double goalanglemax;       /* stop improvement if largest angle reaches this threshold */

		/* quality file output */
		int minsineout;              /* en/disable .minsine file output */
		int minangout;               /* en/disable .minang file output */
		int maxangout;               /* en/disable .maxang file output */
		int vlrmsout;                /* en/disable .vlrms file output */
		int nrrout;                 /* en/disable .rnrr file output */

		/* output file name prefix */
		char fileprefix[100];

		/* enable animation */
		int animate;
		/* for animation, only output .stats */
		int timeseries;

		/* verbosity */
		int verbosity;
		int usecolor;

		/* miscellaneous */
		int outputandquit;           /* just produce all output files for unchanged mesh */
	};

	/* a tet for the improvement stack, referred to by the tuple of its vertices
	* and a single quality measure
	*/
	struct ImproveTetra
	{
		TetraHandle handle;
		double quality;

		bool operator < (ImproveTetra &tet_) const
		{
			if (quality < tet_.quality)
				return true;
			return false;
		}
	};

	typedef std::deque<ImproveTetra> TetraStack;


	/*****************************************************************************/
	/*                                                                           */
	/*  Global variables                                                         */
	/*                                                                           */
	/*****************************************************************************/

	/* quality thresholds for averages */
	extern double meanthresholds[NUMQUALMEASURES][NUMMEANTHRESHOLDS];

	/* types of quality measures that may be used */
	enum TetQualityMetrics
	{
		QUAL_MINSINE,
		QUAL_BIASEDMINSINE,
		QUAL_RADIUSRATIO,
		QUAL_VLRMS3RATIO,
		QUAL_MEANSINE,
		QUAL_MINSINEANDEDGERATIO,
		QUAL_WARPEDMINSINE,
		QUAL_MINANGLE,
		QUAL_MAXANGLE
	};

	enum journalentryclasses
	{
		ALLENTRIES,    /* match all entries */
		INSERTDELETE,  /* insertion or deletion of vertices */
		SMOOTH,        /* change position of a single vertex */
		TOPOLOGICAL,   /* topological change involving multiple vertices */
		LABEL          /* classify a vertex */
	};

	enum journalentrytypes
	{
		INSERTVERTEX,  /* insertion of a vertex */
		DELETEVERTEX,  /* deletion of a vertex */
		SMOOTHVERTEX,  /* smoothing of a vertex */
		DELETETET,     /* single tet deletion */
		INSERTTET,     /* single tet insertion */
		EDGECONTRACT,  /* edge contract */
		FLIP23,        /* 2-3 flip, base case for face removal */
		FLIP32,        /* 3-2 flip, base case for edge removal */
		FLIP22,        /* 2-2 flip */
		FLIP41,        /* 4-1 flip (precedes deletion of interior vertex) */
		FLIP14,        /* 1-4 flip (follows insertion of interior vertex) */
		FLIP13,        /* 1-3 flip (same as 1-4 but inserted vertex is on facet) */
		FLIP31,        /* the inverse of a 1-3 flip */
		FLIP12,        /* 1-2 flip (same as 1-4 but inserted vertex is on edge) */
		FLIP21,        /* inverse of 1-2 */
		CLASSIFY       /* classify vertex (FREEVERTEX, FACETVERTEX, etc) */
	};

	extern double degtosin(double inangle);
}
#endif