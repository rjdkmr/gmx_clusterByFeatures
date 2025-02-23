/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018-2019  Rajendra Kumar
 *
 * gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gmx_clusterByFeatures is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.
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

#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include <iostream>

#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/utility/futil.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/commandline/cmdlineinit.h"

#include "parseData.h"
#include "hole_radius.h"


void CopyRightMsgHole() {
    
    std::string msg = R"~(
        :-)  gmx_clusterByFeatures hole (-:
        
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


int cat_pdb(int nframe, char *fn_input, FILE *fOut)	{
    int i = 0;
    int number=0;
    char **data=NULL;
    
    
    FILE *f_input;
    f_input = fopen(fn_input, "r");
    data = get_all_lines(f_input, &number);
    fclose(f_input);
    
    fprintf(fOut, "\nMODEL    %4d\n", nframe+1);
    for(i=0;i<number;i++)
        fprintf(fOut, "%s", data[i]);
    fprintf(fOut, "TER\n");
    fprintf(fOut, "ENDMDL\n");
    
    free(data);
    return 0;
}

int add_data_to_file(char *fn_input, FILE *fOutResult, rvec cvect)	{
    int i = 0, j = 0, index =0, num;
    FILE *f_input;
    int number=0, dataNumber=0;
    float *radius=NULL, *center=NULL;
    char **data=NULL;
    char **SplitData = NULL;
    gmx_bool bCenXYZ=FALSE;
    gmx_bool bGotDataSection = FALSE;
    float coord[3], **coordList;
    char residues[2][10], **residueList;
    int centerAxis = ZZ, sign=1;
    
    if(cvect[XX] != 0)
        centerAxis = XX;
    
    if(cvect[YY] != 0)
        centerAxis = YY;
    
    if(cvect[ZZ] != 0)
        centerAxis = ZZ;
    
    // when negative vector is used, hole program gives opposite coordinate in radius section
    // so, multiply it by -1
    if(cvect[centerAxis] < 0 )
        sign = -1;
    
    f_input = fopen(fn_input, "r");
    data = get_all_lines(f_input, &number);
    fclose(f_input);
    
    // Parse center coordinate and radius, skip mid-point
    for(i=0;i<number;i++)	{
        
        if (data[i]==NULL)
            continue;
        
        if(strstr(data[i],"cenxyz.cvec      radius")!=NULL)	{
            bCenXYZ = TRUE;
            continue;
        }
        
        if(strstr(data[i],"Minimum radius found")!=NULL)	{
            bCenXYZ = FALSE;
            continue;
        }
        
        if(bCenXYZ)	{
            
            if( (is_first_numerics(data[i])) && (strstr(data[i], "sampled") != NULL) ) {
                
                if(radius == NULL)	{
                    snew(radius, 1);
                    snew(center, 1);
                }
                else	{
                    srenew(radius, dataNumber+1);
                    srenew(center, dataNumber+1);
                }
                
                SplitData = split_by_space(data[i], &num);
                //fprintf(fOutResult,"%s   %s\n",SplitData[0],SplitData[1]);  // old one
                
                dataNumber +=1 ;
                radius[dataNumber-1] = strtof(SplitData[1], NULL);
                center[dataNumber-1] = sign *strtof(SplitData[0], NULL); // multiply by sign
                freeCharArray(SplitData, num);
                SplitData = NULL;
            }
        }
    }
    
    // Allocation memory for residue and coordinate list
    snew(residueList, dataNumber);
    snew(coordList, dataNumber);
    for(i=0; i<dataNumber; i++)	{
        snew(residueList[i], 20);
        snew(coordList[i], 3);
    }
    
    // Parse residues and full coordinate of each center
    for(i=0;i<number;i++)	{
        if (data[i]==NULL)
            continue;
        
        // start for each center
        if(strstr(data[i],"highest radius point found:")!=NULL)	{
            bGotDataSection = TRUE;
            continue;
        }
        
        // finish for each center
        if(strstr(data[i],"stored as")!=NULL)	{
            bGotDataSection = FALSE;
            
            // get residue information for center
            for(j=0; j< dataNumber; j++)	{
                if (coord[centerAxis] == center[j])	{
                    index = j;
                    break;
                }
            }
            
            // Store data in arrays
            sprintf(residueList[index], "%s,%s", residues[0], residues[1]);
            for(j=0; j<3; j++)
                coordList[index][j] = coord[j];
            continue;
        }
        
        
        // Parse coordinate
        if ( (bGotDataSection) && (strstr(data[i],"at point")!=NULL) ) {
            SplitData = split_by_space(data[i], &num);
            coord[XX] = strtof(SplitData[2], NULL);
            coord[YY] = strtof(SplitData[3], NULL);
            coord[ZZ] = strtof(SplitData[4], NULL);
            freeCharArray(SplitData, num);
            SplitData = NULL;
        }
        
        // Parse closest residue
        if ( (bGotDataSection) && (strstr(data[i],"closest atom surface")!=NULL) ) {
            SplitData = split_by_space(data[i], &num);
            if(num == 8)
                sprintf(residues[0], "%s%s", SplitData[5], SplitData[7]);
            else
                sprintf(residues[0], "%s%s", SplitData[5], SplitData[6]);
            freeCharArray(SplitData, num);
            SplitData = NULL;
        }
        
        // Parse 2nd closest residue
        if ( (bGotDataSection) && (strstr(data[i],"2nd closest surface")!=NULL) ) {
            SplitData = split_by_space(data[i], &num);
            // Sometimes ouput contains wiered output, segmentation fault at random
            if(is_first_numerics(SplitData[3])) {
                if(num == 8)
                    sprintf(residues[1], "%s%s", SplitData[5], SplitData[7]);
                else
                    sprintf(residues[1], "%s%s", SplitData[5], SplitData[6]);
            }
            else {
                sprintf(residues[1], "     ");
            }
            freeCharArray(SplitData, num);
            SplitData = NULL;
        }
        
    }
    
    for(i=0; i<dataNumber; i++)	{
        fprintf(fOutResult, "%5.3f  %5.3f  %5.3f  %5.3f  %s\n", coordList[i][0], coordList[i][1], coordList[i][2], radius[i], residueList[i]);
    }
    
    for(i=0; i<dataNumber; i++)	{
        sfree(residueList[i]);
        sfree(coordList[i]);
    }
    sfree(residueList);
    sfree(coordList);
    sfree(radius);
    sfree(center);
    freeCharArray(data, number);
    
    return 0;
}

int calculate_com(t_topology *top, int indsize, int *index, rvec *x, rvec com)	{
    int d, i;
    real mass = 0;
    
    for(i=0;(i<indsize);i++) {
        mass += top->atoms.atom[index[i]].m;
    }
    
    for(d=0;(d<DIM);d++) {
        com[d]=0;
        for(i=0;(i<indsize);i++) {
            com[d] += x[index[i]][d] * top->atoms.atom[index[i]].m * 10;
        }
        com[d] /= mass;
    }
    
    return 0;
}

int gmx_hole (int argc,char *argv[])	{
    const char *desc[] = {
        "To calculate channel radius using hole program. \"hole\" program should",
        "be present in PATH variable. If -fit is enabled,",
        "the molecule will be translated, rotated and centered at the origin. ",
        "First, check with -pdb option for first frame that the hole program ",
        "is able to identify channel correctly. Modify the values of -endrad, ",
        "-sample,-cvect and -cpoint for correct identification of channel with ",
        "-pdb option. The channel could be visualized by using this pdb file "
        "in any visualization program."
    };
    real endrad = 5;
    real sample = 0.5;
    rvec cvect = {0, 0, 1}, cpoint = {999, 999, 999};
    gmx_bool bFit=TRUE;
    int eHoleInputRadius;
    int centatomid=-1;
    int mcstep=1500;
    gmx_output_env_t *oenv;
    
    t_pargs pa[] = {
        { "-fit",    TRUE,  etBOOL, {&bFit},    "To fit structure" },
        { "-endrad", FALSE, etREAL, {&endrad},  "radius above which the hole2 program regards a result as an indicating that the end of the pore has been reached" },
        { "-sample", FALSE, etREAL, {&sample},  "The distance between the planes for the sphere centers." },
        { "-cvect",  FALSE, etRVEC, {cvect},    "This specified a vector which lies in the direction of the channel/pore." },
        { "-cpoint", FALSE, etRVEC, {cpoint},   "A point which lies within the channel and acts as a seed for channel/cavity. If not given, center of mass will be used." },
        { "-catmid",  FALSE, etINT,  {&centatomid},   "Serial number of atom, which lies within the channel and acts as a seed for channel/cavity. If not given, center of mass will be used." },
        { "-mcstep",  FALSE, etINT,  {&mcstep},   "Number of Monte-Carlo steps for hole program" },
        { "-rad",    TRUE,  etENUM, {hole_inp_radius },"Radius of atoms. Accepted categories of radii are" },
    };
    
    t_filenm   fnm[] = {
        { efTRX, "-f",   NULL,      ffREAD },
        { efTPS, NULL,   NULL,      ffREAD },
        { efNDX, NULL,   NULL,      ffOPTRD },
        { efDAT, "-o",  "radius", ffWRITE },
        { efPDB, "-pdb",  "sphpdb", ffOPTWR }
    };
    
    #define NFILE asize(fnm)
    int npargs;
    CopyRightMsgHole();
    
    npargs = asize(pa);
    if ( ! parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT , NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv) )
    {
        return 0;
    }
    
    //Check which radius category is selected 
    eHoleInputRadius = nenum(hole_inp_radius);
    
    FILE *fOutResult, *fOutPDB;
    std::ofstream fHoleRad;
    
    t_trxstatus *status;
    t_topology top;
    //t_atoms    atoms;
    //char title[STRLEN];
    PbcType        ePBC;
    int natoms, nframe=0;
    real       t ;
    matrix     box;
    int 		  indsize, nfit;
    char       *grpnm=NULL,*fitname;
    int    *index=NULL,*ifit=NULL;
    rvec       *xp, *x, x_shift = {0, 0, 0};
    char       prefix_name[32], pdbfile[256], hole_outfile[256], hole_outPDB[256];
    const char *fnOutPDB=NULL;
    char       hole_cmd[1024];
    FILE *tmpf;
    int i=0;
    gmx_bool bOutPDB=FALSE, bM=TRUE;
    
    //Reading topology
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),&top,&ePBC,&xp,NULL,box,FALSE);
    
    if(opt2fn_null("-pdb",NFILE,fnm)!=NULL)	{
        fnOutPDB = opt2fn_null("-pdb",NFILE,fnm);
        bOutPDB = TRUE;
    }
    
    if(bFit)	{
        printf("\nChoose a group for the least squares fit\n");
        get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nfit,&ifit,&fitname);
        if (nfit < 3)
            gmx_fatal(FARGS,"Need >= 3 points to fit!\n");
    }
    
    if (fn2bTPX(ftp2fn(efTPS,NFILE,fnm)) == FALSE)
        bM = FALSE;
    
    real *w_rls=NULL;
    if(bFit)	{
        snew(w_rls,top.atoms.nr);
        for(i=0; (i<nfit); i++)	{
            if(bM)
                w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
            else
                w_rls[ifit[i]]=1.0;
        }
    }
    
    
    //Getting index
    printf("\nChoose a group for hole calculation...\n");
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&indsize,&index,&grpnm);
    
    if (bFit) {
        copy_rvec(xp[index[0]], x_shift);
        reset_x(nfit,ifit,top.atoms.nr,NULL,xp,w_rls);
        rvec_dec(x_shift, xp[index[0]]);
    }
    
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    if (bFit)	{
        reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
        do_fit(natoms,w_rls,xp,x);
        for (i = 0; (i < natoms); i++)
            rvec_inc(x[i], x_shift);
    }
    
    // Check cpoint is provided by the users, if not provided, fill cpoint from center of mass
    if( (cpoint[XX]==999) && (cpoint[YY]==999) && (cpoint[ZZ]==999) )
        calculate_com(&top, indsize, index, x, cpoint);
    
    // Check center-atom-id is given and within the range. If not given, center of mass is already calculated above 
    if(centatomid > 0) {
        gmx_bool bMatchCenterAtomId = FALSE;
        for(i=0; (i<indsize); i++)	{
            if(centatomid-1 == index[i]) {
                bMatchCenterAtomId = TRUE;
            }
        }
        if(!bMatchCenterAtomId) {
            gmx_fatal(FARGS,"Input serial number of atom %d is not found in selected atom-groups index.\n", centatomid);
        }
        printf("\n\nSelected \"%s\" atom of \"%s-%d\" as seed for channel/cavity\n\n", *top.atoms.atomname[centatomid-1], \
        *top.atoms.resinfo[top.atoms.atom[centatomid-1].resind].name, top.atoms.resinfo[top.atoms.atom[centatomid-1].resind].nr);
    }
    
    //Creating name for temporary PDB file
    strcpy(prefix_name,"ddXXXXXX");
    sprintf(prefix_name,"ddXXXXXX");
    gmx_tmpnam(prefix_name);
    remove(prefix_name);
    sprintf(pdbfile,"%s.pdb",prefix_name);
    
    if ((tmpf = fopen(pdbfile,"w")) == NULL)
        gmx_fatal(FARGS,"Can not open pdb file %s",pdbfile);
    fclose(tmpf);
    remove(pdbfile);
    
    //Creating name for output file of hole program
    sprintf(hole_outfile,"%s.out",prefix_name);
    sprintf(hole_outPDB,"%s_sphere.pdb",prefix_name);
    
    // Creating radius file for hole program
    fHoleRad.open("input_atom_radius.rad", std::fstream::out);
    fHoleRad<< holeInputRadiusMap.at(eHoleInputRadius) << std::endl;
    
    
    fOutResult = gmx_ffopen(opt2fn("-o",NFILE,fnm),"w");
    fprintf(fOutResult,"#Axis \t Radius\n");
    
    if(bOutPDB)
        fOutPDB = gmx_ffopen(fnOutPDB,"w");
    
    
    do	{
        if (bFit)	{
            reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
            do_fit(natoms,w_rls,xp,x);
            for (i = 0; (i < natoms); i++)
                rvec_inc(x[i], x_shift);
        }
        
        
        tmpf = gmx_ffopen(pdbfile,"w");
        write_pdbfile_indexed(tmpf,NULL,&top.atoms,x,ePBC,box,' ',-1,indsize,index,NULL,TRUE);
        gmx_ffclose(tmpf);
        
        // Assign cpoint from input atom coordinate
        if(centatomid > 0) {
            cpoint[XX] = x[centatomid-1][XX]*10;
            cpoint[YY] = x[centatomid-1][YY]*10;
            cpoint[ZZ] = x[centatomid-1][ZZ]*10;
        }

        //Creating variable for executing hole
        if(bOutPDB)
            sprintf(hole_cmd,"hole >%s <<EOF\ncoord %s\nradius %s\ncvect %2.2f %2.2f %2.2f\nsample %2.2f\nendrad %2.2f\ncpoint %2.2f %2.2f %2.2f\nmcstep %d\nsphpdb %s\nEOF", \
            hole_outfile, pdbfile, "input_atom_radius.rad", cvect[XX], cvect[YY], cvect[ZZ], sample, endrad, cpoint[XX], cpoint[YY], cpoint[ZZ], mcstep, hole_outPDB );
        else
            sprintf(hole_cmd,"hole >%s <<EOF\ncoord %s\nradius %s\ncvect %2.2f %2.2f %2.2f\nsample %2.2f\nendrad %2.2f\ncpoint %2.2f %2.2f %2.2f\n\nmcstep %d\nEOF", \
            hole_outfile, pdbfile,  "input_atom_radius.rad", cvect[XX], cvect[YY], cvect[ZZ], sample, endrad, cpoint[XX], cpoint[YY], cpoint[ZZ], mcstep);
        
        if(0 != system(hole_cmd))
            gmx_fatal(FARGS,"Failed to execute command: %s",hole_cmd);
        
        fprintf(fOutResult,"\n# Time = %15.5f\n",t);
        add_data_to_file(hole_outfile, fOutResult, cvect);
        
        if(bOutPDB)
            cat_pdb(nframe, hole_outPDB , fOutPDB);
        
        remove(pdbfile);
        remove(hole_outPDB);
        nframe++;
        
    }while(read_next_x(oenv,status,&t,x,box));
    
    
    fprintf(stdout, "Thanks for using gmx_hole!!!\n");
    fprintf(stdout, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
    fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
    fprintf(stderr, "O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)\n");
    fprintf(stderr, "The Pore Dimensions of Gramicidin A\n");
    fprintf(stderr, "Biophysical Journal 65:2455-2460.\n");
    fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
    
    return 0;
}

int main(int argc, char *argv[])
{
    
    #ifdef GMX_NO_SYSTEM
    gmx_fatal(FARGS,"No calls to system(3) supported on this platform.");
    #endif
    
    gmx_run_cmain(argc, argv, &gmx_hole);
    
    return 0;
}
