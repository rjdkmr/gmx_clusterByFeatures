/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018-2025  Rajendra Kumar
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


#ifndef HOLE_RADIUS_H
#define HOLE_RADIUS_H

#include<string>
#include <map>

    
std::string amberuni = R"~(remark: Input file for program Tooshort.
remark: Contains bond and vdw radius records for
remark: for pdb atoms.
remark: format:
remark:   BOND C??? 0.8 
remark:   = bond radius of atom with first character C
remark:   is 0.8 angs.  ? is a wildcard character
remark: keyword VDWR for van der Waals radii
remark:    VDWR CA   GLY 1.925
remark: N.B. PUT SPECIFIC RECORDS BEFORE GENERAL
remark:
remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
remark: van der Waals radii
remark: all cb's are type C2 except thr, val, ile
VDWR CB   ALA 2.00
VDWR CB   THR 1.85
VDWR CB   VAL 1.85
VDWR CB   ILE 1.85
VDWR CB   ??? 1.925
remark: other C2 atoms (carbon with two aliphatic hydrogens)
VDWR CA   GLY 1.925
VDWR CG   GLU 1.925
VDWR CG   LYS 1.925
VDWR CD   LYS 1.925
VDWR CE   LYS 1.925
VDWR CG   PRO 1.925
VDWR CD   PRO 1.925
VDWR CG   MET 1.925
VDWR CG1  ILE 1.925
VDWR CG   GLN 1.925
VDWR CG   ARG 1.925
VDWR CD   ARG 1.925
remark: C3 atoms (carbon with two aliphatic hydrogens)
VDWR CD?  LEU 2.00
VDWR CE   MET 2.00
VDWR CG2  THR 2.00
VDWR CD1  ILE 2.00
VDWR CG2  ILE 2.00
VDWR CG?  VAL 2.00
remark: NH3 atom type N3
VDWR NZ   LYS 1.85
remark: sp2 oxygen atom type O
VDWR OE1  GLN 1.60
VDWR O    ??? 1.60
remark: acid oxygens O2
VDWR OE?  GLU 1.60
VDWR OD?  ASP 1.60
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR P??? ??? 2.10)~";

std::string bondi = R"~(remark: van der Waals radii: from Bondi, A. (1964) J. Phys. Chem. 68:441-451 
remark: shorter than amber - take your pick
VDWR C??? ??? 1.70
VDWR O??? ??? 1.52
VDWR S??? ??? 1.80
VDWR N??? ??? 1.55
VDWR H??? ??? 1.20
VDWR P??? ??? 1.80)~";

std::string downscaled = R"~(remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.00
VDWR O??? ??? 0.891891892
VDWR S??? ??? 1.081081081
VDWR N??? ??? 0.945945946
VDWR H??? ??? 0.540540541
VDWR P??? ??? 1.135135135
VDWR Z??? ??? 0.4
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00 
VDWR D2?  ASN 1.00 
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85)~";
 
std::string hardcore = R"~(remark "hard core" radii given by Turano, Pear & Busath (1992)
remark  biophysical journal, vol 63, 152-161 
VDWR C??? ??? 1.425
VDWR O??? ??? 1.404
VDWR N??? ??? 1.28
VDWR H??? ??? 0.576
remark rest scaled from simple.rad by 1.425/1.85=0.77
VDWR S??? ??? 1.54
VDWR P??? ??? 1.617
VDWR ???? ??? 1.4
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85)~";
 
std::string simple = R"~(remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR P??? ??? 2.10
VDWR Z??? ??? 0.74
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00 
VDWR D2?  ASN 1.00 
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85)~";

std::string xplor = R"~(remark: van der Waals radii from XPLOR parameter file:param19x.pro
remark: VDWR determined as sigma/2
remark: Kindly provided by Ian Kerr, Molecular Biophysics, Oxford University
remark:
remark: CB ALA = CH3E, CB THR, VAL, ILE = CH1E
remark: CB ??? = CH2E
VDWR CB   ALA 1.929
VDWR CB   THR 2.107
VDWR CB   VAL 2.107
VDWR CB   ILE 2.107
VDWR CB   ??? 1.991
remark: other CH2E atoms (carbon with two aliphatic hydrogens)
VDWR CA   GLY 1.991
VDWR CG   GLU 1.991
VDWR CG   LYS 1.991
VDWR CD   LYS 1.991
VDWR CE   LYS 1.991
VDWR CG   PRO 1.991
VDWR CD   PRO 1.991
VDWR CG   MET 1.991
VDWR CG1  ILE 1.991
VDWR CG   GLN 1.991
VDWR CG   ARG 1.991
VDWR CD   ARG 1.991
remark: CH3E atoms (carbon with three aliphatic hydrogens)
remark: smaller than CH2E ? I'm only copying from Xplor
VDWR CD?  LEU 1.929
VDWR CE   MET 1.929
VDWR CG2  THR 1.929
VDWR CD1  ILE 1.929
VDWR CG2  ILE 1.929
VDWR CG?  VAL 1.929
remark: aromatic carbon atoms
VDWR CE3  TRP 1.871
VDWR CD1  TRP 1.871
VDWR CZ?  TRP 1.871
VDWR CH2  TRP 1.871
VDWR CD?  PHE 1.871
VDWR CE?  PHE 1.871
VDWR CD?  TYR 1.871
VDWR CE?  TYR 1.871
VDWR CD2  HIS 1.871
VDWR CE1  HIS 1.871
VDWR CZ   ??? 1.871
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 0.713
VDWR D2?  ASN 0.713
remark: lysine terminal hydrogens
VDWR HZ?  LYS 0.535
remark: general last, all oxygens have same vdwr in XPLOR
remark: all sulphurs have same vdwr in XPLOR
remark: all nitrogens same vdwr IN XPLOR
remark: all other H will be polar - same vdwr
remark: remaining C, eg. CA, carbonyl C
VDWR C??? ??? 1.871
VDWR H??? ??? 0.713
VDWR O??? ??? 1.425
VDWR S??? ??? 1.684
VDWR N??? ??? 1.425
VDWR P??? ??? 1.800)~";

const char *hole_inp_radius[] = { NULL,  "bondi", "amberuni", "downscaled", "hardcore", "simple", "xplor", NULL };
enum { eBondiHoleRadius =1, eAmberuniHoleRadius, eDownscaledHoleRadius, eHardcoreHoleRadius, eSimpleHoleRadius, eXplorHoleRadius };
std::map<int, std::string> holeInputRadiusMap = {
    {eAmberuniHoleRadius, amberuni},
    {eBondiHoleRadius, bondi},
    {eDownscaledHoleRadius, downscaled},
    {eHardcoreHoleRadius, hardcore},
    {eSimpleHoleRadius, simple},
    {eXplorHoleRadius, xplor},    
};


#endif

