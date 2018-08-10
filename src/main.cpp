
#include "gromacs/commandline/cmdlineinit.h"

int gmx_clusterByFeatures(int argc,char *argv[]);

int main (int argc,char *argv[])	{
    gmx_run_cmain(argc,argv, &gmx_clusterByFeatures);
	return 0;
}
