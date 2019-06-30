/*
 * This file is part of g_distMat
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014-2015  Rajendra Kumar
 *
 * g_distMat is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_distMat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_distMat.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <cmath>
#include <iostream>
#include <cstring>
#include <pthread.h>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>

#include "gromacs/commandline/viewit.h"
#include "gromacs/utility/futil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/index.h"
#include "gromacs/pbcutil/rmpbc.h"

std::vector<std::string> split(const std::string &s, char delim);

void *status;
pthread_t *thread;
pthread_attr_t attr;
pthread_mutex_t dist_mutex;
int ithread;

typedef struct {
	int nA, nB, nframes, nPcaCoords;
	int *resndxA, *resndxB, *natmresA, *natmresB;
	real cutoff, **mean, **var, **cmap, **dist, **sumdist, **sumsqdist;
	rvec *coord, *pcaCoords;
	int nthreads;
	gmx_bool b2ndTraj;
    matrix box = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
} DIST_MAT;

DIST_MAT distance_matrix;

void CopyRightMsgDistMat() {

    std::string msg = R"~(
         :-)  gmx_clusterByFeatures distmat (-:

             Author: Rajendra Kumar

       Copyright (C) 2018-2019  Rajendra Kumar


gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gmx_clusterByFeatures is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    )~";

        std::cerr<<msg<<"\n";
}



void set_num_threads(int nt)	{
	distance_matrix.nthreads = nt;
}


void make_thread_index()	{
	long i = 0, start, end;

	int NTHREADS = distance_matrix.nthreads;

	for(i=0; i<NTHREADS; i++){
		if ((distance_matrix.nA%NTHREADS) < 1)	{
			start = ceil(i*((distance_matrix.nA/NTHREADS)));
			end = start + ceil(distance_matrix.nA/NTHREADS);
		}
		else	{
			start = ceil(i*((distance_matrix.nA/NTHREADS)));
			end = start + ceil(distance_matrix.nA/NTHREADS);
			if (i==NTHREADS-1)
				end = start + ceil(distance_matrix.nA/NTHREADS)+distance_matrix.nA%NTHREADS;
		}
		//end = ceil((i+1)*(distance_matrix.nA/NTHREADS))+NTHREADS;
		fprintf(stdout,"\n%ld %ld\n",start, end);
	}
}

// Function to calculate minimum distance between two residues
real calc_min_dist (int resi, int *resndxA, int *natmresA, int resj, int *resndxB, int *natmresB, rvec *x)
{
	real dist=999.0, r;
	int i, j,ii,jj;


	for(i=0;i<natmresA[resi];i++)
	{
		ii = resndxA[resi] + i;
		for(j=0;j<natmresB[resj];j++)
		{
			jj = resndxB[resj] + j;

			r = ((x[ii][XX]-x[jj][XX])*(x[ii][XX]-x[jj][XX]))+
			    ((x[ii][YY]-x[jj][YY])*(x[ii][YY]-x[jj][YY])) +
			    ((x[ii][ZZ]-x[jj][ZZ])*(x[ii][ZZ]-x[jj][ZZ]));

			r = sqrt(r);

			if (r<dist)
				dist = r;
		}
	}

	return dist;
}

// Initialization of residue indexing for group A
void init_ndx_distmat_grp_A(t_topology top, int *index, int isize)	{
	int i=0, nres=0, prevres;

	prevres = -1;
	nres  = 0;

	snew(distance_matrix.natmresA, 1);
	snew(distance_matrix.resndxA, 1);

	for(i=0;i<isize;i++)
	{
		if(top.atoms.atom[index[i]].resind == prevres)	{
			distance_matrix.natmresA[nres-1]++;
		}

		if(top.atoms.atom[index[i]].resind != prevres)
		{
			nres++;

			srenew(distance_matrix.natmresA, nres);
			srenew(distance_matrix.resndxA, nres);

			distance_matrix.natmresA[nres-1] = 0;
			distance_matrix.natmresA[nres-1]++;
			distance_matrix.resndxA[nres-1] = index[i];

			prevres = top.atoms.atom[index[i]].resind;
		}
	}
	distance_matrix.nA = nres;
	fprintf(stderr,"There are %d residues with %d atoms in first group\n",nres,isize);
}

// Initialization of residue indexing for group B
void init_ndx_distmat_grp_B(t_topology top, int *index, int isize)	{
	int i=0, nres=0, prevres;

	prevres = -1;
	nres  = 0;

	snew(distance_matrix.natmresB, 1);
	snew(distance_matrix.resndxB, 1);

	for(i=0;i<isize;i++)
	{
		if(top.atoms.atom[index[i]].resind == prevres)	{
			distance_matrix.natmresB[nres-1]++;
		}

		if(top.atoms.atom[index[i]].resind != prevres)
		{
			nres++;

			srenew(distance_matrix.natmresB, nres);
			srenew(distance_matrix.resndxB, nres);

			distance_matrix.natmresB[nres-1] = 0;
			distance_matrix.natmresB[nres-1]++;
			distance_matrix.resndxB[nres-1] = index[i];

			prevres = top.atoms.atom[index[i]].resind;
		}
	}
	distance_matrix.nB = nres;
	fprintf(stderr,"There are %d residues with %d atoms in second group\n",nres,isize);
}

// Initialization of output data
void init_outdata_distmat(int gapX, int gapY, bool index_grp_same)	{
	int i=0, j=0, d=0, m=0, nB=0;

    snew(distance_matrix.dist, distance_matrix.nA);
	snew(distance_matrix.sumdist, distance_matrix.nA);
	snew(distance_matrix.sumsqdist, distance_matrix.nA);
	snew(distance_matrix.mean, distance_matrix.nA);
	snew(distance_matrix.var, distance_matrix.nA);
	snew(distance_matrix.cmap, distance_matrix.nA);


	for(i=0; i<distance_matrix.nA; i++){
        snew(distance_matrix.dist[i], distance_matrix.nB);
		snew(distance_matrix.sumdist[i], distance_matrix.nB);
		snew(distance_matrix.sumsqdist[i], distance_matrix.nB);
		snew(distance_matrix.mean[i], distance_matrix.nB);
		snew(distance_matrix.var[i], distance_matrix.nB);
		snew(distance_matrix.cmap[i], distance_matrix.nB);


		for(j=0; j<distance_matrix.nB; j++){
            distance_matrix.dist[i][j] = 0;
			distance_matrix.sumdist[i][j] = 0;
			distance_matrix.sumsqdist[i][j] = 0;
			distance_matrix.mean[i][j] = 0;
			distance_matrix.var[i][j] = 0;
			distance_matrix.cmap[i][j] = 0;
		}

	}
	distance_matrix.nframes = 0;
    
    // Initialize PCA-coords vector
    i = 0;
    j = 0;
    d = 0;
    m = 1;
    distance_matrix.nPcaCoords = 1;
    while (i< distance_matrix.nA) {
        
        if (index_grp_same)
            nB = i;
        else
            nB = distance_matrix.nB;
        
        while (j < nB) {
            d++; // XX, YY, or ZZ
            m++; // add total count of elements
            if (d == DIM)
            {
                d = 0;
                distance_matrix.nPcaCoords += 1 ;
            }
            j = j + gapY;
        }
        i = i + gapX;
        j = 0 ;
    }
    snew(distance_matrix.pcaCoords, distance_matrix.nPcaCoords);
    fprintf(stdout, " Number of distance-matrix elements for PCA trajectory: %d\n", m);
    fprintf(stdout, " Number of distance-matrix coordinates in PCA trajectory: %d\n", distance_matrix.nPcaCoords);
}

void calculate_dist_mat() {
	int i=0, j=0;
	real dist=0;

	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			 dist = calc_min_dist(i, distance_matrix.resndxA, distance_matrix.natmresA, j, distance_matrix.resndxB, distance_matrix.natmresB, distance_matrix.coord);
             distance_matrix.dist[i][j] = dist;
			 distance_matrix.sumdist[i][j] += dist;
			 distance_matrix.sumsqdist[i][j] += (dist*dist);
             if (dist <= distance_matrix.cutoff) {
                 distance_matrix.cmap[i][j]++;
             }
		}
	}
}


void *calculate_dist_mat_pthread (void *arg) {
	int NTHREADS = distance_matrix.nthreads;
	long t;
	t = (long)arg;
	int start, end;

	if ((distance_matrix.nA%NTHREADS) < 1)	{
		start = ceil(t*((distance_matrix.nA/NTHREADS)));
		end = start + ceil(distance_matrix.nA/NTHREADS);
	}
	else	{
		start = ceil(t*((distance_matrix.nA/NTHREADS)));
		end = start + ceil(distance_matrix.nA/NTHREADS);
		if (t==NTHREADS-1)
			end = start + ceil(distance_matrix.nA/NTHREADS)+distance_matrix.nA%NTHREADS;
	}

	int i=start, j=0;
	real dist=0;

	if(NTHREADS==1){
		start = 0;
		end = distance_matrix.nA;
	}
	if(end>distance_matrix.nA)
		end=distance_matrix.nA;


	for(i=start; i<end; i++){
		for(j=0; j<distance_matrix.nB; j++){
			 dist = calc_min_dist(i, distance_matrix.resndxA, distance_matrix.natmresA, j, distance_matrix.resndxB, distance_matrix.natmresB, distance_matrix.coord);
			 pthread_mutex_lock(&dist_mutex);

			 // For first trajectory
			 if (!distance_matrix.b2ndTraj)		{
                 distance_matrix.dist[i][j] = dist;
				 distance_matrix.sumdist[i][j] += dist;
				 distance_matrix.sumsqdist[i][j] += (dist*dist);

				 if (dist <= distance_matrix.cutoff)
					 distance_matrix.cmap[i][j]++;
			 }
			 else	{
                 distance_matrix.dist[i][j] = dist;
				 distance_matrix.sumdist[i][j] += dist;
				 distance_matrix.sumsqdist[i][j] += ( (distance_matrix.mean[i][j] - dist) * (distance_matrix.mean[i][j] - dist) );
			 }

			 pthread_mutex_unlock(&dist_mutex);
		}
	}
	if(NTHREADS>1)
		pthread_exit(NULL);
	return (void *)0;
}

void calc_mean_var()	{
	int i=0, j=0, n = distance_matrix.nframes;
	real sum=0, sumsq=0;
	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			sum = distance_matrix.sumdist[i][j];
			sumsq = distance_matrix.sumsqdist[i][j];

			distance_matrix.mean[i][j] = sum/n ;
			distance_matrix.var[i][j] = (sumsq - ( (sum*sum) /n)) / (n-1);
			distance_matrix.cmap[i][j] /= n;

		}

	}
}

void calc_msd_2nd()	{
	int i=0, j=0, n = distance_matrix.nframes;
	real sum=0, sumsq=0;
	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			sum = distance_matrix.sumdist[i][j];
			sumsq = distance_matrix.sumsqdist[i][j];

			distance_matrix.mean[i][j] = abs( (sum/n) - distance_matrix.mean[i][j] );
			distance_matrix.var[i][j] = sumsq / n;

		}

	}
}

void write_distmat_xtc(struct t_fileio *fioXTC, int gapX, int gapY, int step, real time, bool index_grp_same, real power) {
    
    int              i=0, j=0, d=0, l=0, nB=0, ret_val;
    real            prec=1000;
    
    while (i< distance_matrix.nA) {
        
        if (index_grp_same)
            nB = i;
        else
            nB = distance_matrix.nB;

        while (j < nB) {
            if (power > 0) {
                if (power == 1)    {
                    distance_matrix.pcaCoords[l][d] = distance_matrix.dist[i][j];
                }
                else    {
                    distance_matrix.pcaCoords[l][d] = std::pow(distance_matrix.dist[i][j], power);
                }
            }
            else if (power < 0){
                if (power == -1)    {
                    distance_matrix.pcaCoords[l][d] = 1/distance_matrix.dist[i][j];
                }
                else {
                    distance_matrix.pcaCoords[l][d] = 1/std::pow(distance_matrix.dist[i][j], power*-1);
                }
            } else{
                distance_matrix.pcaCoords[l][d] = distance_matrix.dist[i][j];
            }
            d++; // XX, YY, or ZZ
            if (d == DIM)
            {
                d = 0;
                l += 1 ;
            }
            j = j + gapY;
        }
        i = i + gapX;
        j = 0;
    }
    
    ret_val = write_xtc(fioXTC, distance_matrix.nPcaCoords, step, time, distance_matrix.box, distance_matrix.pcaCoords, prec);
    if(!ret_val) {
        gmx_fatal(FARGS, "Cannot write trajectory file for PCA.");
    }
}


bool is_index_group_same(int *indexA, int isizeA, int *indexB, int isizeB) {
    int i=0;
    bool same = true;
    if(isizeA == isizeB) {
        for(i=0; i < isizeA; i++) {
            if(indexA[i] != indexB[i]) {
                same = false;
                break;
            }
        }
    }
    else {
        same = false;
    }
    return same;
}

void write_pca_dummy_pdb(t_topology top, int *indexA, int isizeA, std::string fnPcaPdb, rvec *x, int ePBC, matrix box){
    int *outindex, i= 0, j=0;
    FILE *fPDB = gmx_ffopen(fnPcaPdb.c_str(), "w");
    
    snew(outindex, distance_matrix.nPcaCoords);
    while( i < distance_matrix.nPcaCoords){
        for(j = 0; j < isizeA; j++ ){
            outindex[i] = indexA[j];
            i++;
            if (i == distance_matrix.nPcaCoords)
                break;
        }
    }
    
    write_pdbfile_indexed(fPDB, "Dummy PDB for dist-mat PCA", &top.atoms,
                           x, ePBC, box, 'A',
                           1, distance_matrix.nPcaCoords, outindex,
                           NULL, FALSE);
}

int gmx_distMat(int argc,char *argv[])
{
	const char *desc[] = {
			"It calculates average minimum-distance matrix of residues between two atom-groups.",
			"Simultaneously, it calculates variance  and standard-deviation matrices.",
			"Also, it calculates contact-frequency map over the trajectory for the ",
			"residues that are within a minimum distance given by \"-ct\" option value.\n\n",
			"To speed up the calculation, it uses all available cores of the CPU using",
			"multi-threading. Number of threads/cores could be change by \"-nt\" option.\n",
			"\n To calculate distance variances or deviation (fluctuations) in a trajectory",
			" with respect to average distances from another trajectory, use [-f traj_for_average.xtc] ",
			" and [-f2 traj_for_variance.xtc]. The averages will be calculated from ",
			"\"traj_for_average.xtc\". Subsequently, variances and deviation will be ",
			"calculated for \"traj_for_variance.xtc\" with respect to previosly calculated averages.\n",
	};
	static real cutoff = 0.4;
	real power = 1;
	int NTHREADS = sysconf( _SC_NPROCESSORS_ONLN );
    int gapX = 5, gapY = 1;
	t_pargs pa[] = {
			{ "-ct",   FALSE, etREAL, {&cutoff}, "cut-off distance (nm) for contact map" },
			{ "-nt",   FALSE, etINT, {&NTHREADS}, "number of threads for multi-threading" },
            { "-gx",   FALSE, etINT, {&gapX}, "Gap between residues along X-axis in distance-matrix for PCA trajectory" },
            { "-gy",   FALSE, etINT, {&gapY}, "ap between residues along Y-axis in distance-matrix for PCA trajectory" },
            { "-power", FALSE, etREAL, {&power}, "Distances will be raised by this power and then dumped in xtc file."} 	};

	t_filenm   fnm[] = {
			{ efTRX, "-f",  NULL, ffREAD },
			{ efTPS, NULL,  NULL, ffREAD },
			{ efNDX, NULL,  NULL, ffOPTRD },
			{ efTRX, "-f2",  NULL, ffOPTRD },
			{ efDAT, "-mean", "average", ffWRITE },
			{ efDAT, "-var", "variance", ffOPTWR },
			{ efDAT, "-std", "stdeviation", ffOPTWR },
			{ efDAT, "-cmap", "contact_map", ffOPTWR },
			{ efXTC, "-pca", "pca", ffOPTWR }
	};

	#define NFILE asize(fnm)

	FILE *fMean, *fVar=NULL, *fStd=NULL, *fCmap=NULL;
	t_topology top;
	int ePBC;

	int isizeA, isizeB;
	int *indexA, *indexB;
	char *grpnameA, *grpnameB;
    bool index_grp_same;

	int i, j, trxnat, nframe=0;
	long nt;
	t_trxstatus *trjstatus;
	rvec *x;
	char title[256];
    
    struct t_fileio *fioXTC=NULL;
    std::string fnPcaTraj, fnPcaPdb;

	real time;
	matrix box;
	gmx_output_env_t *oenv;
	gmx_rmpbc_t gpbc = NULL;

	// Copyright message
	CopyRightMsgDistMat();

	// Parse command line argument and print all options
		if ( ! parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv) )	{
			return 0;
		}

	// Set number of threads
	set_num_threads(NTHREADS);

	// Set cut-off distance for contact map
	distance_matrix.cutoff = cutoff;

	// Reading tpr filr
	read_tps_conf(ftp2fn(efTPS,NFILE,fnm), &top, &ePBC, &x, NULL, box, FALSE);

	// Selection of first index group
	fprintf(stderr,"Select first group:\n");
	get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isizeA,&indexA,&grpnameA);

	// Selection of second index group
	fprintf(stderr,"Select second group:\n");
	get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isizeB,&indexB,&grpnameB);
    
    // Check if index group is same
    index_grp_same = is_index_group_same(indexA, isizeA, indexB, isizeB);

	// Initialization for first group
	init_ndx_distmat_grp_A(top, indexA, isizeA);

	// Initialization for second group
	init_ndx_distmat_grp_B(top, indexB, isizeB);

	// Initialization for output arrays
	init_outdata_distmat(gapX, gapY, index_grp_same);
    
    // Open xtc file
    if (opt2bSet("-pca", NFILE, fnm)) {
        std::vector < std::string > temp;
        fnPcaTraj =  opt2fn_null("-pca", NFILE, fnm);
        fioXTC = open_xtc(opt2fn("-pca", NFILE, fnm), "w");
        temp = split(fnPcaTraj, '.');
        fnPcaPdb = temp[temp.size()-2] + "_dummy.pdb";
    }

	trxnat=read_first_x(oenv, &trjstatus,ftp2fn(efTRX,NFILE,fnm), &time, &distance_matrix.coord, box);

	gpbc = gmx_rmpbc_init(&top.idef, ePBC, trxnat);

	thread = (pthread_t *) malloc(sizeof(pthread_t) * NTHREADS);

	do {
		distance_matrix.nframes++;

		if(NTHREADS>1)	{

			//make_thread_index();
			//exit(1);

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			pthread_mutex_init(&dist_mutex, NULL);

			for(nt=0;nt<NTHREADS;nt++)
				ithread = pthread_create(&thread[nt],&attr, calculate_dist_mat_pthread, (void *) nt);
				//mi_calculate((void *)&i);

			 //Free attribute and wait for the other threads
			pthread_attr_destroy(&attr);
			for(nt=0; nt<NTHREADS; nt++) {
				ithread = pthread_join(thread[nt], &status);
				if (ithread) {
					printf("ERROR; return code from pthread_join() is %d\n", ithread);
					exit(-1);
				 }
			 }
		}

		else
            calculate_dist_mat();
        
        if (fioXTC != NULL) {
            write_distmat_xtc(fioXTC, gapX, gapY, nframe, time, index_grp_same, power);
        }
        
        nframe++;
        
	}while (read_next_x(oenv,trjstatus, &time, distance_matrix.coord, box));

    // Write pdb file for PCA
    if (fioXTC != NULL) {
        write_pca_dummy_pdb(top, indexA, isizeA, fnPcaPdb, x, ePBC, box);
    }
    
	// Final calculation of average, var, std-deviation and contact-map
	calc_mean_var();

	// For 2nd trajectory /////////////////////////////////
	if (opt2fn_null("-f2", NFILE, fnm) != NULL )	{
		distance_matrix.b2ndTraj = TRUE;
	}

	if (distance_matrix.b2ndTraj)	{
		fprintf(stdout, "\nReading second trajectory ...\n");

		for(i=0; i<distance_matrix.nA; i++)	{
			for(j=0; j<distance_matrix.nB; j++)	{
				distance_matrix.sumdist[i][j] = 0;
				distance_matrix.sumsqdist[i][j] = 0;
			}
		}

		trxnat=read_first_x(oenv, &trjstatus,opt2fn_null("-f2", NFILE, fnm), &time, &distance_matrix.coord, box);
		gpbc = gmx_rmpbc_init(&top.idef, ePBC, trxnat);

		do {

			if(NTHREADS>1)	{

				//make_thread_index();
				//exit(1);

				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				pthread_mutex_init(&dist_mutex, NULL);

				for(nt=0;nt<NTHREADS;nt++)
					ithread = pthread_create(&thread[nt],&attr, calculate_dist_mat_pthread, (void *) nt);
					//mi_calculate((void *)&i);

				 //Free attribute and wait for the other threads
				pthread_attr_destroy(&attr);
				for(nt=0; nt<NTHREADS; nt++) {
					ithread = pthread_join(thread[nt], &status);
					if (ithread) {
						printf("ERROR; return code from pthread_join() is %d\n", ithread);
						exit(-1);
					 }
				 }
			}

			else
				calculate_dist_mat();
            
		} while (read_next_x(oenv,trjstatus, &time, distance_matrix.coord, box));
	}
	//////////////////////////////////////////////////////////


	//Writing output file
	fMean = opt2FILE("-mean", NFILE, fnm, "w");
    
    if (opt2bSet("-var", NFILE, fnm))    {
        fVar = opt2FILE("-var", NFILE, fnm, "w");
    }
    
    if (opt2bSet("-std", NFILE, fnm))    {
        fStd = opt2FILE("-std", NFILE, fnm, "w");
    }

	if ((!distance_matrix.b2ndTraj) && (opt2bSet("-cmap", NFILE, fnm)) )    {
        fCmap = opt2FILE("-cmap", NFILE, fnm, "w");
    }


	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			fprintf(fMean, "%5.3f ", distance_matrix.mean[i][j]);
            
            if (fVar != NULL)   fprintf(fVar, "%5.3f ", distance_matrix.var[i][j]);
            
            if (fStd != NULL)   fprintf(fStd, "%5.3f ", sqrt(distance_matrix.var[i][j]));

			if ( (!distance_matrix.b2ndTraj) && (fCmap != NULL) )
				fprintf(fCmap, "%5.12f ", distance_matrix.cmap[i][j]);
		}
		fprintf(fMean,"\n");
		if (fVar != NULL) fprintf(fVar,"\n");
		if (fStd != NULL) fprintf(fStd,"\n");

		if ( (!distance_matrix.b2ndTraj) && (fCmap != NULL) )
			fprintf(fCmap,"\n");
	}
	
	return 0;
}
