#define MAXINCIDENTPOINT 150
#define MAXINCIDENTHETRDRON 300
#define HUGEFLOAT 1.0e20
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define MAXVERTEX 0xFFFF
//define quality kind
#define MINSINE 0
#define VOLLENGTH 1
#define RADIUSRATIO 2
#define MINSINE2 3
#define BIASEDMINSINE 4
#define MAXCHANGE 100

#define ACTIVESETFACTOR 1.03
#define MAXACTIVESET 150
#define MAXBASIS 100
#define NEARESTMIN 1.0e-13
/* rate must be worse by this much to reset alpha */
#define RATEEPSILON 1.0e-6
/* minimum quality improvement in a smoothing step */
#define MINSMOOTHITERIMPROVE 1.0e-5
/* if d is within this of zero-length, call it zero */
#define DEPSILON 1.0e-5
/* minimum step size */
#define MINSTEPSIZE 1.0e-5
/* maximum iterations in non-smooth line search */
#define MAXLINEITER 50
/* maximum iterations of non smooth optimization */
#define MAXSMOOTHITER 5
#define SMOOTHTIME 1

//for cuda
#define CUDAHUGEFLOAT 1000000000
#define CUDANEARESTMIN 0.0000000000001
/* rate must be worse by this much to reset alpha */
#define CUDARATEEPSILON 0.000001
/* minimum quality improvement in a smoothing step */
#define CUDAMINSMOOTHITERIMPROVE 0.00001
/* if d is within this of zero-length, call it zero */
#define CUDADEPSILON 0.00001
/* minimum step size */
#define CUDAMINSTEPSIZE 0.00001