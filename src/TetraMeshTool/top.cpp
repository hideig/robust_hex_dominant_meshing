#include <TetraMeshTool/top.h>
namespace VolumeMesh
{
	/* quality thresholds for averages */
	double meanthresholds[NUMQUALMEASURES][NUMMEANTHRESHOLDS] =
	{
		/* QUALMINSINE thresholds */
		{SINE1,SINE5,SINE10,SINE15,SINE25,SINE35,SINE45},
		/* QUALRADIUSRATIO thresholds */
		{0.1,0.2,0.3,0.4,0.5,0.6,0.7},
		/* QUALVLRMS3RATIO thresholds */
		{0.1,0.2,0.3,0.4,0.5,0.6,0.7},
		/* QUALMEANSINE thresholds */
		{SINE1,SINE5,SINE10,SINE15,SINE25,SINE35,SINE45},
		/* QUALMINSINEANDEDGERATIO thresholds */
		{SINE1,SINE5,SINE10,SINE15,SINE25,SINE35,SINE45},
		/* QUALWARPEDMINSINE thresholds */
		{SINE1,SINE5,SINE10,SINE15,SINE25,SINE35,SINE45}
	};


	/* convert from an angle in degrees to sin */
	double degtosin(double inangle)
	{
		return sin((inangle) * (PI / 180.0));
	}


}