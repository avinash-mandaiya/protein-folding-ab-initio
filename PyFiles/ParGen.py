#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:49:34 2023

@author: avinashmandaiya
"""
import os
import numpy as np    
import sys
from Bio.PDB import *
import _pickle as pickle
from numpy import random
import matplotlib.pyplot as plt
import math
from Bio.PDB import vectors
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.vectors import calc_angle
from Bio.PDB.vectors import Vector
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.DSSP import DSSP
import warnings
from sklearn.decomposition import PCA
import timeit
# import editdistance

warnings.resetwarnings()
warnings.filterwarnings(action="ignore", message="unclosed", category=ResourceWarning,append=True)
warnings.filterwarnings(action="error",append=True)

#%%

def residue_type(AA):
    AA_list = "RNDQEHKCMGPSTWYILVFA"
    result = AA_list.find(AA)
    
    return result

def trans_origin(ini_data):
    no_rows = ini_data.shape[0]
    avg = ini_data.sum(axis=0)/no_rows
    avg_matrix = np.tile(avg, (no_rows,1)) 
    ini_data = ini_data - avg_matrix
    
    return ini_data

#%%

lenAA = np.zeros(20,dtype = int)
    
lenAA[0] = 11
lenAA[1] = 8  
lenAA[2] = 8  
lenAA[3] = 9  
lenAA[4] = 9  
lenAA[5] = 10 
lenAA[6] = 9  
lenAA[7] = 6  
lenAA[8] = 8  
lenAA[9] = 5  
lenAA[10] = 7  
lenAA[11] = 6  
lenAA[12] = 7  
lenAA[13] = 14 
lenAA[14] = 12 
lenAA[15] = 8  
lenAA[16] = 8  
lenAA[17] = 7  
lenAA[18] = 11 
lenAA[19] = 5 

AA_list = "RNDQEHKCMGPSTWYILVFA"

#%%
"""rigid motifs (side-chain and peptide bond)"""

"""fragments of side chains"""

"""disulphide bonds: SS = 0 means no disulphide bond"""
SS = 1

numF = np.zeros(20,dtype = int)

numF[0] = 5
numF[1] = 3
numF[2] = 3
numF[3] = 4
numF[4] = 4
numF[5] = 3
numF[6] = 6

if SS == 0:
    numF[7] = 3
else:
    numF[7] = 2
    
numF[8] = 5
numF[9] = 1
numF[10] = 2
numF[11] = 3
numF[12] = 4
numF[13] = 3
numF[14] = 3
numF[15] = 5
numF[16] = 5
numF[17] = 4
numF[18] = 3
numF[19] = 2

sizeF = np.zeros((21,np.max(numF)),dtype = int)

sizeF[0][0:numF[0]] = [5,5,5,5,10]
sizeF[1][0:numF[1]] = [5,5,6]
sizeF[2][0:numF[2]] = [5,5,4]
sizeF[3][0:numF[3]] = [5,5,5,6]
sizeF[4][0:numF[4]] = [5,5,5,4]
# sizeF[5][0:numF[5]] = [5,5,10]
sizeF[5][0:numF[5]] = [5,5,9]
sizeF[6][0:numF[6]] = [5,5,5,5,5,5]

if SS == 0:
    sizeF[7][0:numF[7]] = [5,5,3]
else:
    sizeF[7][0:numF[7]] = [5,5]
    
sizeF[8][0:numF[8]] = [5,5,5,3,5]
sizeF[9][0:numF[9]] = [5]
sizeF[10][0:numF[10]] = [5,11]
sizeF[11][0:numF[11]] = [5,5,3]
sizeF[12][0:numF[12]] = [5,5,3,5]
sizeF[13][0:numF[13]] = [5,5,16]
sizeF[14][0:numF[14]] = [5,5,13]
sizeF[15][0:numF[15]] = [5,5,5,5,5]
sizeF[16][0:numF[16]] = [5,5,5,5,5]
sizeF[17][0:numF[17]] = [5,5,5,5]
sizeF[18][0:numF[18]] = [5,5,12]
sizeF[19][0:numF[19]] = [5,5]


atom_numF = np.zeros((20,np.max(numF),np.max(sizeF)),dtype = int)

atom_numF[0][0][0:sizeF[0][0]] = [0,1,2,4,11]
atom_numF[0][1][0:sizeF[0][1]] = [1,4,5,12,13]
atom_numF[0][2][0:sizeF[0][2]] = [4,5,6,14,15]
atom_numF[0][3][0:sizeF[0][3]] = [5,6,7,16,17]
atom_numF[0][4][0:sizeF[0][4]] = [6,7,8,9,10,18,19,20,21,22]

atom_numF[1][0][0:sizeF[1][0]] = [0,1,2,4,8]
atom_numF[1][1][0:sizeF[1][1]] = [1,4,5,9,10]
atom_numF[1][2][0:sizeF[1][2]] = [4,5,6,7,11,12]

atom_numF[2][0][0:sizeF[2][0]] = [0,1,2,4,8]
atom_numF[2][1][0:sizeF[2][1]] = [1,4,5,9,10]
atom_numF[2][2][0:sizeF[2][2]] = [4,5,6,7]

atom_numF[3][0][0:sizeF[3][0]] = [0,1,2,4,9]
atom_numF[3][1][0:sizeF[3][1]] = [1,4,5,10,11]
atom_numF[3][2][0:sizeF[3][2]] = [4,5,6,12,13]
atom_numF[3][3][0:sizeF[3][3]] = [5,6,7,8,14,15]

atom_numF[4][0][0:sizeF[4][0]] = [0,1,2,4,9]
atom_numF[4][1][0:sizeF[4][1]] = [1,4,5,10,11]
atom_numF[4][2][0:sizeF[4][2]] = [4,5,6,12,13]
atom_numF[4][3][0:sizeF[4][3]] = [5,6,7,8]

atom_numF[5][0][0:sizeF[5][0]] = [0,1,2,4,10]
atom_numF[5][1][0:sizeF[5][1]] = [1,4,5,11,12]
# atom_numF[5][2][0:sizeF[5][2]] = [4,5,6,7,8,9,13,14,15,16]
atom_numF[5][2][0:sizeF[5][2]] = [4,5,6,7,8,9,13,14,15]

atom_numF[6][0][0:sizeF[6][0]] = [0,1,2,4,9]
atom_numF[6][1][0:sizeF[6][1]] = [1,4,5,10,11]
atom_numF[6][2][0:sizeF[6][2]] = [4,5,6,12,13]
atom_numF[6][3][0:sizeF[6][3]] = [5,6,7,14,15]
atom_numF[6][4][0:sizeF[6][4]] = [6,7,8,16,17]
atom_numF[6][5][0:sizeF[6][5]] = [7,8,18,19,20]

if SS == 0:
    atom_numF[7][0][0:sizeF[7][0]] = [0,1,2,4,6]
    atom_numF[7][1][0:sizeF[7][1]] = [1,4,5,7,8]
    atom_numF[7][2][0:sizeF[7][2]] = [4,5,9]
else:
    atom_numF[7][0][0:sizeF[7][0]] = [0,1,2,4,6]
    atom_numF[7][1][0:sizeF[7][1]] = [1,4,5,7,8]

atom_numF[8][0][0:sizeF[8][0]] = [0,1,2,4,8]
atom_numF[8][1][0:sizeF[8][1]] = [1,4,5,9,10]
atom_numF[8][2][0:sizeF[8][2]] = [4,5,6,11,12]
atom_numF[8][3][0:sizeF[8][3]] = [5,6,7]
atom_numF[8][4][0:sizeF[8][4]] = [6,7,13,14,15]

atom_numF[9][0][0:sizeF[9][0]] = [0,1,2,4,5]

atom_numF[10][0][0:sizeF[10][0]] = [0,1,2,4,7]
atom_numF[10][1][0:sizeF[10][1]] = [0,1,4,5,6,8,9,10,11,12,13]

atom_numF[11][0][0:sizeF[11][0]] = [0,1,2,4,6]
atom_numF[11][1][0:sizeF[11][1]] = [1,4,5,7,8]
atom_numF[11][2][0:sizeF[11][2]] = [4,5,9]

atom_numF[12][0][0:sizeF[12][0]] = [0,1,2,4,7]
atom_numF[12][1][0:sizeF[12][1]] = [1,4,5,6,8]
atom_numF[12][2][0:sizeF[12][2]] = [4,5,9]
atom_numF[12][3][0:sizeF[12][3]] = [4,6,10,11,12]

atom_numF[13][0][0:sizeF[13][0]] = [0,1,2,4,14]
atom_numF[13][1][0:sizeF[13][1]] = [1,4,5,15,16]
atom_numF[13][2][0:sizeF[13][2]] = [4,5,6,7,8,9,10,11,12,13,17,18,19,20,21,22]

atom_numF[14][0][0:sizeF[14][0]] = [0,1,2,4,12]
atom_numF[14][1][0:sizeF[14][1]] = [1,4,5,13,14]
atom_numF[14][2][0:sizeF[14][2]] = [4,5,6,7,8,9,10,11,15,16,17,18,19]

atom_numF[15][0][0:sizeF[15][0]] = [0,1,2,4,8]
atom_numF[15][1][0:sizeF[15][1]] = [1,4,5,6,9]
atom_numF[15][2][0:sizeF[15][2]] = [4,5,7,10,11]
atom_numF[15][3][0:sizeF[15][3]] = [4,6,12,13,14]
atom_numF[15][4][0:sizeF[15][4]] = [5,7,15,16,17]

atom_numF[16][0][0:sizeF[16][0]] = [0,1,2,4,8]
atom_numF[16][1][0:sizeF[16][1]] = [1,4,5,9,10]
atom_numF[16][2][0:sizeF[16][2]] = [4,5,6,7,11]
atom_numF[16][3][0:sizeF[16][3]] = [5,6,12,13,14]
atom_numF[16][4][0:sizeF[16][4]] = [5,7,15,16,17]

atom_numF[17][0][0:sizeF[17][0]] = [0,1,2,4,7]
atom_numF[17][1][0:sizeF[17][1]] = [1,4,5,6,8]
atom_numF[17][2][0:sizeF[17][2]] = [4,5,9,10,11]
atom_numF[17][3][0:sizeF[17][3]] = [4,6,12,13,14]

atom_numF[18][0][0:sizeF[18][0]] = [0,1,2,4,11]
atom_numF[18][1][0:sizeF[18][1]] = [1,4,5,12,13]
atom_numF[18][2][0:sizeF[18][2]] = [4,5,6,7,8,9,10,14,15,16,17,18]

atom_numF[19][0][0:sizeF[19][0]] = [0,1,2,4,5]
atom_numF[19][1][0:sizeF[19][1]] = [1,4,6,7,8]

#%%
RAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CG"},
    "CD" : {"charge" : 0.20, "charmmid" : "CT2", "atom_name" : "CD"},
    "NE" : {"charge" : -0.70, "charmmid" : "NC2", "atom_name" : "NE"},
    "CZ" : {"charge" : 0.64, "charmmid" : "C", "atom_name" : "CZ"},
    "NH1" : {"charge" : -0.80, "charmmid" : "NC2", "atom_name" : "NH1"},
    "NH2" : {"charge" : -0.80, "charmmid" : "NC2", "atom_name" : "NH2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "HD2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD1"},
    "HD3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD2"},
    "HE" : {"charge" : 0.44, "charmmid" : "HC", "atom_name" : "HE"},
    "HH11" : {"charge" : 0.46, "charmmid" : "HC", "atom_name" : "HH11"},
    "HH12" : {"charge" : 0.46, "charmmid" : "HC", "atom_name" : "HH12"},
    "HH21" : {"charge" : 0.46, "charmmid" : "HC", "atom_name" : "HH21"},
    "HH22" : {"charge" : 0.46, "charmmid" : "HC", "atom_name" : "HH22"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

NAA = {  
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},       
    "CG" : {"charge" : 0.55, "charmmid" : "CC", "atom_name" : "CG"},
    "OD1" : {"charge" : -0.55, "charmmid" : "O", "atom_name" : "OD1"},
    "ND2" : {"charge" : -0.62, "charmmid" : "NH2", "atom_name" : "ND2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HD21" : {"charge" : 0.32, "charmmid" : "H", "atom_name" : "HD21"},
    "HD22" : {"charge" : 0.30, "charmmid" : "H", "atom_name" : "HD22"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

DAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.28, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : 0.62, "charmmid" : "CC", "atom_name" : "CG"},
    "OD1" : {"charge" : -0.76, "charmmid" : "OC", "atom_name" : "OD1"},
    "OD2" : {"charge" : -0.76, "charmmid" : "OC", "atom_name" : "OD2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

QAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CG"},
    "CD" : {"charge" : 0.55, "charmmid" : "CC", "atom_name" : "CD"},
    "OE1" : {"charge" : -0.55, "charmmid" : "O", "atom_name" : "OE1"},
    "NE2" : {"charge" : -0.62, "charmmid" : "NH2", "atom_name" : "NE2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "HE21" : {"charge" : 0.32, "charmmid" : "H", "atom_name" : "HE21"},
    "HE22" : {"charge" : 0.30, "charmmid" : "H", "atom_name" : "HE22"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

EAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.28, "charmmid" : "CT2", "atom_name" : "CG"},
    "CD" : {"charge" : 0.62, "charmmid" : "CC", "atom_name" : "CD"},
    "OE1" : {"charge" : -0.76, "charmmid" : "OC", "atom_name" : "OE1"},
    "OE2" : {"charge" : -0.76, "charmmid" : "OC", "atom_name" : "OE2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

# HAA = {    
#     "N" : {"charge" : -0.47, "charmmid" : "NH1"},
#     "CA" : {"charge" : 0.07, "charmmid" : "CT1"},
#     "C" : {"charge" : 0.51, "charmmid" : "C"},
#     "O" : {"charge" : -0.51, "charmmid" : "O"},
#     "CB" : {"charge" : -0.05, "charmmid" : "CT2"},
#     "CG" : {"charge" : 0.19, "charmmid" : "CPH1"},
#     "ND1" : {"charge" : -0.51, "charmmid" : "NR3"},
#     "CD2" : {"charge" : 0.19, "charmmid" : "CPH1"},
#     "CE1" : {"charge" : 0.32, "charmmid" : "CPH2"},
#     "NE2" : {"charge" : -0.51, "charmmid" : "NR3"},
#     "HA" : {"charge" : 0.09, "charmmid" : "HB1"},
#     "HB2" : {"charge" : 0.09, "charmmid" : "HA"},
#     "HB3" : {"charge" : 0.09, "charmmid" : "HA"},
#     "HD1" : {"charge" : 0.44, "charmmid" : "H"},
#     "HD2" : {"charge" : 0.13, "charmmid" : "HR1"},
#     "HE1" : {"charge" : 0.18, "charmmid" : "HR2"},
#     "HE2" : {"charge" : 0.44, "charmmid" : "H"},
#     "H" : {"charge" : 0.31, "charmmid" : "H"}
#     }

HAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.09, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.05, "charmmid" : "CPH1", "atom_name" : "CG"},
    "ND1" : {"charge" : -0.36, "charmmid" : "NR1", "atom_name" : "ND1"},
    "CD2" : {"charge" : 0.22, "charmmid" : "CPH1", "atom_name" : "CD2"},
    "CE1" : {"charge" : 0.25, "charmmid" : "CPH2", "atom_name" : "CE1"},
    "NE2" : {"charge" : -0.70, "charmmid" : "NR2", "atom_name" : "NE2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HD1" : {"charge" : 0.32, "charmmid" : "H", "atom_name" : "HD1"},
    "HD2" : {"charge" : 0.10, "charmmid" : "HR3", "atom_name" : "HD2"},
    "HE1" : {"charge" : 0.13, "charmmid" : "HR1", "atom_name" : "HE1"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

KAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CG"},
    "CD" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CD"},
    "CE" : {"charge" : 0.21, "charmmid" : "CT2", "atom_name" : "CE"},
    "NZ" : {"charge" : -0.30, "charmmid" : "NH3", "atom_name" : "NZ"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "HD2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD1"},
    "HD3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD2"},
    "HE2" : {"charge" : 0.05, "charmmid" : "HA2", "atom_name" : "HE1"},
    "HE3" : {"charge" : 0.05, "charmmid" : "HA2", "atom_name" : "HE2"},
    "HZ1" : {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HZ1"},
    "HZ2" : {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HZ2"},
    "HZ3" : {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HZ3"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

if SS == 0:
    CAA = {    
        "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
        "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
        "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
        "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
        "CB" : {"charge" : -0.11, "charmmid" : "CT2", "atom_name" : "CB"},
        "SG" : {"charge" : -0.23, "charmmid" : "S", "atom_name" : "SG"},
        "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
        "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
        "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
        "HG" : {"charge" : 0.16, "charmmid" : "HS", "atom_name" : "HG"},
        "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
        }
else:
    CAA = {    
        "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
        "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
        "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
        "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
        "CB" : {"charge" : -0.10, "charmmid" : "CT2", "atom_name" : "CB"},
        "SG" : {"charge" : -0.08, "charmmid" : "S", "atom_name" : "SG"},
        "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
        "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
        "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
        "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
        }    

MAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.14, "charmmid" : "CT2", "atom_name" : "CG"},
    "SD" : {"charge" : -0.09, "charmmid" : "S", "atom_name" : "SD"},
    "CE" : {"charge" : -0.22, "charmmid" : "CT3", "atom_name" : "CE"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "HE1" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HE1"},
    "HE2" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HE2"},
    "HE3" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HE3"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

GAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : -0.02, "charmmid" : "CT2", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "HA2" : {"charge" : 0.09, "charmmid" : "HB2", "atom_name" : "HA1"},
    "HA3" : {"charge" : 0.09, "charmmid" : "HB2", "atom_name" : "HA2"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

PAA = {    
    "N" : {"charge" : -0.29, "charmmid" : "N", "atom_name" : "N"},
    "CA" : {"charge" : 0.02, "charmmid" : "CP1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CP2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.18, "charmmid" : "CP2", "atom_name" : "CG"},
    "CD" : {"charge" : 0.00, "charmmid" : "CP3", "atom_name" : "CD"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG1"},
    "HG3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG2"},
    "HD2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD1"},
    "HD3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HD2"}
    }

SAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : 0.05, "charmmid" : "CT2", "atom_name" : "CB"},
    "OG" : {"charge" : -0.66, "charmmid" : "OH1", "atom_name" : "OG"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG" : {"charge" : 0.43, "charmmid" : "H", "atom_name" : "HG1"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

TAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : 0.14, "charmmid" : "CT1", "atom_name" : "CB"},
    "OG1" : {"charge" : -0.66, "charmmid" : "OH1", "atom_name" : "OG1"},
    "CG2" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CG2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB" : {"charge" : 0.09, "charmmid" : "HA1", "atom_name" : "HB"},
    "HG1" : {"charge" : 0.43, "charmmid" : "H", "atom_name" : "HG1"},
    "HG21" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG21"},
    "HG22" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG22"},
    "HG23" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG23"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

WAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.03, "charmmid" : "CY", "atom_name" : "CG"},
    "CD1" : {"charge" : 0.035, "charmmid" : "CA", "atom_name" : "CD1"},
    "CD2" : {"charge" : -0.02, "charmmid" : "CPT", "atom_name" : "CD2"},
    "NE1" : {"charge" : -0.61, "charmmid" : "NY", "atom_name" : "NE1"},
    "CE2" : {"charge" : 0.13, "charmmid" : "CPT", "atom_name" : "CE2"},
    "CE3" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CE3"},
    "CZ2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CZ2"},
    "CZ3" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CZ3"},
    "CH2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CH2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HD1" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HD1"},
    "HE1" : {"charge" : 0.38, "charmmid" : "H", "atom_name" : "HE1"},
    "HE3" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HE3"},
    "HZ2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HZ2"},
    "HZ3" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HZ3"},
    "HH2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HH2"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

YAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : 0.00, "charmmid" : "CA", "atom_name" : "CG"},
    "CD1" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CD1"},
    "CD2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CD2"},
    "CE1" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CE1"},
    "CE2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CE2"},
    "CZ" : {"charge" : 0.11, "charmmid" : "CA", "atom_name" : "CZ"},
    "OH" : {"charge" : -0.54, "charmmid" : "OH1", "atom_name" : "OH"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HD1" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HD1"},
    "HD2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HD2"},
    "HE1" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HE1"},
    "HE2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HE2"},
    "HH" : {"charge" : 0.43, "charmmid" : "H", "atom_name" : "HH"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

IAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.09, "charmmid" : "CT1", "atom_name" : "CB"},
    "CG1" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CG1"},
    "CG2" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CG2"},
    "CD1" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CD"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB" : {"charge" : 0.09, "charmmid" : "HA1", "atom_name" : "HB"},
    "HG12" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG11"},
    "HG13" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HG12"},
    "HG21" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG21"},
    "HG22" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG22"},
    "HG23" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG23"},
    "HD11" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD1"},
    "HD12" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD2"},
    "HD13" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD3"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

LAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : -0.09, "charmmid" : "CT1", "atom_name" : "CG"},
    "CD1" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CD1"},
    "CD2" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CD2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HG" : {"charge" : 0.09, "charmmid" : "HA1", "atom_name" : "HG"},
    "HD11" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD11"},
    "HD12" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD12"},
    "HD13" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD13"},
    "HD21" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD21"},
    "HD22" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD22"},
    "HD23" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HD23"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

VAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.09, "charmmid" : "CT1", "atom_name" : "CB"},
    "CG1" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CG1"},
    "CG2" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CG2"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB" : {"charge" : 0.09, "charmmid" : "HA1", "atom_name" : "HB"},
    "HG11" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG11"},
    "HG12" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG12"},
    "HG13" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG13"},
    "HG21" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG21"},
    "HG22" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG22"},
    "HG23" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HG23"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

FAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.18, "charmmid" : "CT2", "atom_name" : "CB"},
    "CG" : {"charge" : 0.00, "charmmid" : "CA", "atom_name" : "CG"},
    "CD1" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CD1"},
    "CD2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CD2"},
    "CE1" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CE1"},
    "CE2" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CE2"},
    "CZ" : {"charge" : -0.115, "charmmid" : "CA", "atom_name" : "CZ"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB1"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA2", "atom_name" : "HB2"},
    "HD1" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HD1"},
    "HD2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HD2"},
    "HE1" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HE1"},
    "HE2" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HE2"},
    "HZ" : {"charge" : 0.115, "charmmid" : "HP", "atom_name" : "HZ"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

AAA = {    
    "N" : {"charge" : -0.47, "charmmid" : "NH1", "atom_name" : "N"},
    "CA" : {"charge" : 0.07, "charmmid" : "CT1", "atom_name" : "CA"},
    "C" : {"charge" : 0.51, "charmmid" : "C", "atom_name" : "C"},
    "O" : {"charge" : -0.51, "charmmid" : "O", "atom_name" : "O"},
    "CB" : {"charge" : -0.27, "charmmid" : "CT3", "atom_name" : "CB"},
    "HA" : {"charge" : 0.09, "charmmid" : "HB1", "atom_name" : "HA"},
    "HB1" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HB1"},
    "HB2" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HB2"},
    "HB3" : {"charge" : 0.09, "charmmid" : "HA3", "atom_name" : "HB3"},
    "H" : {"charge" : 0.31, "charmmid" : "H", "atom_name" : "HN"}
    }

#%%

res_info = {
    "R" : dict(RAA),
    "N" : dict(NAA),
    "D" : dict(DAA),
    "Q" : dict(QAA),
    "E" : dict(EAA),
    "H" : dict(HAA),
    "K" : dict(KAA),
    "C" : dict(CAA),
    "M" : dict(MAA),
    "G" : dict(GAA),
    "P" : dict(PAA),
    "S" : dict(SAA),
    "T" : dict(TAA),
    "W" : dict(WAA),
    "Y" : dict(YAA),
    "I" : dict(IAA),
    "L" : dict(LAA),
    "V" : dict(VAA),
    "F" : dict(FAA),
    "A" : dict(AAA) 
    }

atom_num_AA = np.zeros(20, dtype= int)
for i in range(20):
    aa = AA_list[i]
    atom_num_AA[i] = len(res_info[aa])

#%%
res_info_all = {
    "R" : dict(RAA),
    "N" : dict(NAA),
    "D" : dict(DAA),
    "Q" : dict(QAA),
    "E" : dict(EAA),
    "H" : dict(HAA),
    "K" : dict(KAA),
    "C" : dict(CAA),
    "M" : dict(MAA),
    "G" : dict(GAA),
    "P" : dict(PAA),
    "S" : dict(SAA),
    "T" : dict(TAA),
    "W" : dict(WAA),
    "Y" : dict(YAA),
    "I" : dict(IAA),
    "L" : dict(LAA),
    "V" : dict(VAA),
    "F" : dict(FAA),
    "A" : dict(AAA), 
 
    "RNT" : dict(RAA),
    "NNT" : dict(NAA),
    "DNT" : dict(DAA),
    "QNT" : dict(QAA),
    "ENT" : dict(EAA),
    "HNT" : dict(HAA),
    "KNT" : dict(KAA),
    "CNT" : dict(CAA),
    "MNT" : dict(MAA),
    "GNT" : dict(GAA),
    "PNT" : dict(PAA),
    "SNT" : dict(SAA),
    "TNT" : dict(TAA),
    "WNT" : dict(WAA),
    "YNT" : dict(YAA),
    "INT" : dict(IAA),
    "LNT" : dict(LAA),
    "VNT" : dict(VAA),
    "FNT" : dict(FAA),
    "ANT" : dict(AAA), 

    "RCT" : dict(RAA),
    "NCT" : dict(NAA),
    "DCT" : dict(DAA),
    "QCT" : dict(QAA),
    "ECT" : dict(EAA),
    "HCT" : dict(HAA),
    "KCT" : dict(KAA),
    "CCT" : dict(CAA),
    "MCT" : dict(MAA),
    "GCT" : dict(GAA),
    "PCT" : dict(PAA),
    "SCT" : dict(SAA),
    "TCT" : dict(TAA),
    "WCT" : dict(WAA),
    "YCT" : dict(YAA),
    "ICT" : dict(IAA),
    "LCT" : dict(LAA),
    "VCT" : dict(VAA),
    "FCT" : dict(FAA),
    "ACT" : dict(AAA)
    }

#%%
for key in res_info_all.keys():
    
    if len(key) == 3:
        
        if key[-2:] == 'NT':
            if key[0] == 'P':
                (res_info_all[key])['H1'] = {"charge" : 0.24, "charmmid" : "HC", "atom_name" : "HN1"}
                (res_info_all[key])['H2'] = {"charge" : 0.24, "charmmid" : "HC", "atom_name" : "HN2"}
                res_info_all[key].update({"N" : {"charge" : -0.07, "charmmid" : "NP", "atom_name" : "N"}})
                res_info_all[key].update({"CA" : {"charge" : 0.16, "charmmid" : "CP1", "atom_name" : "CA"}})
                res_info_all[key].update({"CD" : {"charge" : 0.16, "charmmid" : "CP3", "atom_name" : "CD"}})

            elif key[0] == 'G':
                (res_info_all[key])['H1'] = (res_info_all[key]).pop('H')
                res_info_all[key].update({"H1" : {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT1"}})
                (res_info_all[key])['H2'] = {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT2"}
                (res_info_all[key])['H3'] = {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT3"}
                res_info_all[key].update({"N" : {"charge" : -0.30, "charmmid" : "NH3", "atom_name" : "N"}})
                res_info_all[key].update({"CA" : {"charge" : 0.13, "charmmid" : "CT2", "atom_name" : "CA"}})

            else:
                (res_info_all[key])['H1'] = (res_info_all[key]).pop('H')
                res_info_all[key].update({"H1" : {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT1"}})
                (res_info_all[key])['H2'] = {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT2"}
                (res_info_all[key])['H3'] = {"charge" : 0.33, "charmmid" : "HC", "atom_name" : "HT3"}
                res_info_all[key].update({"N" : {"charge" : -0.30, "charmmid" : "NH3", "atom_name" : "N"}})
                res_info_all[key].update({"CA" : {"charge" : 0.21, "charmmid" : "CT1", "atom_name" : "CA"}})
                res_info_all[key].update({"HA" : {"charge" : 0.10, "charmmid" : "HB1", "atom_name" : "HA"}})
 
        if key[-2:] == 'CT':
            res_info_all[key]['OXT'] = {"charge" : -0.67, "charmmid" : "OC", "atom_name" : "OT2"}            
            res_info_all[key].update({"C" : {"charge" : 0.34, "charmmid" : "CC", "atom_name" : "C"}})
            res_info_all[key].update({"O" : {"charge" : -0.67, "charmmid" : "OC", "atom_name" : "OT1"}})
            
#%%
"""get van der waal parameters data"""

full_path = os.path.realpath(__file__)
os.chdir(os.path.dirname(full_path))
os.chdir('../GenFiles')

file1 = 'VDW.par'
vdw_dict = {}

with open(file1, 'r') as f:
    text = f.read()
    
list_text = text.split("\n")

for i in range(len(list_text) - 1):
    x = list_text[i].split()
    if len(x) == 4:
        vdw_dict[x[0]] = {"eps" : float(x[2]), "rmin" : float(x[3])}
    elif len(x) == 7:
        vdw_dict[x[0]] = {"eps" : float(x[2]), "rmin" : float(x[3]), "eps14" : float(x[5]), "rmin14" : float(x[6])}
    # print(len(x))
  
#%%
"""Only SEQ file is available"""

os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')   

seq_name = "1LB0"

filename = seq_name+".seq" 

with open(filename, 'r') as fh:
    Lines = fh.readlines()
    
    pro_length = int((Lines[0].split())[0])
    given_seq = (Lines[1].split())[0]
 

pq = np.zeros((pro_length,26),dtype = float)
atom_yn = np.zeros((pro_length,26),dtype = int)
atom_num_pro = np.zeros(pro_length,dtype = int)

for i in range(pro_length):
    aa = given_seq[i]
    
    atom_num_pro[i] = atom_num_AA[residue_type(aa)]
    
    if i == 0:
        aa += "NT"
        atom_num_pro[i] += 2
    if i == pro_length - 1:
        aa += "CT"
        atom_num_pro[i] += 1
        
    for atom_type in list(res_info_all[aa].keys()):
        j = list(res_info_all[aa].keys()).index(atom_type)

        pq[i,j] = float(res_info_all[aa][atom_type]["charge"])
        atom_yn[i][j] = 1        

#%%
"""PDB file is available"""

""""""""""""""""""""""""
seq_name = "1LB0"
""""""""""""""""""""""""

"""main-protein"""
os.chdir(os.path.dirname(full_path))
os.chdir('../PDB')

filename = seq_name+".pdb"     

mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)

parser = PDBParser()

structure = parser.get_structure(seq_name, filename)  

model = structure[0]
# dssp = DSSP(model, filename, dssp='mkdssp')    # for the dssp module

id_chain = list(model.child_dict.keys())[0]
main_chain = model[id_chain]


all_atoms = Selection.unfold_entities(structure, "A")
all_residues = Selection.unfold_entities(main_chain, "R")

ppb=PPBuilder()

given_seq = ""
for pp in ppb.build_peptides(structure):
    given_seq += str(pp.get_sequence())
    
pro_length = 0

for residue in all_residues:
    if (residue.get_parent()).get_id() == id_chain and residue.get_id()[0] == ' ':
        pro_length = pro_length+1

atom_coor = np.zeros((pro_length,26,3),dtype = float)
pq = np.zeros((pro_length,26),dtype = float)
atom_yn = np.zeros((pro_length,26),dtype = int)
atom_num_pro = np.zeros(pro_length,dtype = int)

for i in range(pro_length):
    res = all_residues[i]
    if (res.get_parent()).get_id() != id_chain:
        continue

    if res.get_id()[0] != ' ':
        continue  
    
    a_key = (id_chain,res.get_id())
    aa = given_seq[i]
    
    atom_num_pro[i] = atom_num_AA[residue_type(aa)]
    
    if i == 0:
        aa += "NT"
        atom_num_pro[i] += 2
    if i == pro_length - 1:
        aa += "CT"
        atom_num_pro[i] += 1
    
    for atom_type in list(res_info_all[aa].keys()):
        j = list(res_info_all[aa].keys()).index(atom_type)
        
        if atom_type[0] == "H":
            try:
                atom_coor[i,j,:] = res[atom_type].get_coord()
            except:
                atom_type_new = atom_type.replace("H", "D", 1)
                atom_coor[i,j,:] = res[atom_type_new].get_coord()
                
        else:
            atom_coor[i,j,:] = res[atom_type].get_coord()
        
        # if atom_type[0] == "H":
        #     continue        

        pq[i,j] = float(res_info_all[aa][atom_type]["charge"])
        if i == pro_length-1:
            print(aa,atom_type,pq[i,j])
        atom_yn[i][j] = 1
                
#%%        

"""read bond data"""

os.chdir(os.path.dirname(full_path))
os.chdir('../GenFiles')
file1 = "AA.bonds"

BondsAA = np.zeros((20,26,26),dtype = float)

with open(file1, 'r') as f:
    text = f.read()
    
list_text = text.split("\n\n")

for i in range(len(list_text)):
    x = list_text[i].split("\n")
    y = x[0].split()
    
    if SS == 0 and y[0] == 'CS':
        continue
    
    if SS == 1 and y[0] == 'CX':
        continue
    
    if len(y[0]) == 2:
        AA = residue_type(y[0][0])
    elif len(y[0]) == 1:    
        AA = residue_type(y[0])
        
    num_bonds = int(y[1])
    
    for j in range(1,num_bonds+1):
        y = x[j].split()
        p = int(y[0])
        q = int(y[1])
        BondsAA[AA,p,q] = 1
        BondsAA[AA,q,p] = 1
 
#%%
tot_atoms = 0

for i in range(pro_length):
    tot_atoms += atom_num_pro[i]
    print(atom_num_pro[i])
    
atom2res_map = np.zeros((tot_atoms,2),dtype = int)
res2atom_map = np.zeros((pro_length,26),dtype = int)

count = 0
for i in range(pro_length):
    for m in range(atom_num_pro[i]):
        atom2res_map[count][0] = i
        atom2res_map[count][1] = m
        res2atom_map[i][m] = count
        count += 1
    

#%%
bond_graph = np.zeros((tot_atoms,tot_atoms),dtype = int)
int_yn = -1*np.ones((pro_length,pro_length,26,26),dtype = int)

for i in range(pro_length):
    for j in range(pro_length):
        for m in range(atom_num_pro[i]):
            for n in range(atom_num_pro[j]):
                if m == n and i == j:
                    continue
                else:
                    int_yn[i][j][m][n] = 0
                    int_yn[j][i][n][m] = 0

"""checking if two atoms belong to the same rigid motif in a particular residue"""

# for i in range(pro_length):
#     AA = residue_type(given_seq[i])
#     for m in range(atom_num_pro[i]):
#         for n in range(m+1,atom_num_pro[i]):
#             for k in range(numF[AA]):
#                 if m in atom_numF[AA][k][0:sizeF[AA][k]] and n in atom_numF[AA][k][0:sizeF[AA][k]]:
#                     int_yn[i][i][m][n] = 4
#                     int_yn[i][i][n][m] = 4

# """cis-trans planar rigid motif"""

# for i in range(pro_length-1):
#     AA2 = residue_type(given_seq[i+1])
    
#     if AA2 == 10:
#         listNH = [0,1]
#     else:
#         listNH = [0,1,atom_num_AA[AA2]-1]
        
#     listCO = [1,2,3]
    
#     for m in range(atom_num_pro[i]):
#         for n in range(atom_num_pro[i+1]):
#             if m in listCO and n in listNH:
#                 int_yn[i][i+1][m][n] = 4
#                 int_yn[i+1][i][n][m] = 4
                
#     for m in range(atom_num_pro[i+1]):
#         for n in range(atom_num_pro[i+1]):
#             if m in listNH and n in listNH:
#                 int_yn[i+1][i+1][m][n] = 4
#                 int_yn[i+1][i+1][n][m] = 4   
                
#     for m in range(atom_num_pro[i]):
#         for n in range(atom_num_pro[i]):                
#             if m in listCO and n in listCO:
#                 int_yn[i][i][m][n] = 4
#                 int_yn[i][i][n][m] = 4  
                
# AA2 = residue_type(given_seq[0])                
# listNH = [0,1,atom_num_pro[0]-3,atom_num_pro[0]-2,atom_num_pro[0]-1]  

# for m in range(atom_num_pro[0]):
#     for n in range(atom_num_pro[0]):
#         if m in listNH and n in listNH:
#             int_yn[0][0][m][n] = 4
#             int_yn[0][0][n][m] = 4   

# AA2 = residue_type(given_seq[pro_length-1])                
# listCO = [1,2,3,atom_num_pro[pro_length-1]-1]       

# for m in range(atom_num_pro[pro_length-1]):
#     for n in range(atom_num_pro[pro_length-1]):                
#         if m in listCO and n in listCO:
#             int_yn[pro_length-1][pro_length-1][m][n] = 4
#             int_yn[pro_length-1][pro_length-1][n][m] = 4   
            
"""Covalent bonds"""

for i in range(pro_length):
    for j in range(i,pro_length):
        if j == i:
            for m in range(atom_num_pro[i]):
                for n in range(m+1,atom_num_pro[j]):
                    AA = residue_type(given_seq[i])
                    
                    if BondsAA[AA,m,n] == 1:
                        int_yn[i][j][m][n] = 1
                        int_yn[j][i][n][m] = 1
                        bond_graph[res2atom_map[i][m]][res2atom_map[j][n]] = 1
                        bond_graph[res2atom_map[j][n]][res2atom_map[i][m]] = 1
             
            if i == 0:
                m = 0
                n = atom_num_pro[i]-2
                int_yn[i][j][m][n] = 1
                int_yn[j][i][n][m] = 1
                bond_graph[res2atom_map[i][m]][res2atom_map[j][n]] = 1
                bond_graph[res2atom_map[j][n]][res2atom_map[i][m]] = 1
                
                m = 0
                n = atom_num_pro[i]-1
                int_yn[i][j][m][n] = 1
                int_yn[j][i][n][m] = 1
                bond_graph[res2atom_map[i][m]][res2atom_map[j][n]] = 1
                bond_graph[res2atom_map[j][n]][res2atom_map[i][m]] = 1
                
            if i == pro_length-1:
                m = 2
                n = atom_num_pro[i]-1
                int_yn[i][j][m][n] = 1
                int_yn[j][i][n][m] = 1
                bond_graph[res2atom_map[i][m]][res2atom_map[j][n]] = 1
                bond_graph[res2atom_map[j][n]][res2atom_map[i][m]] = 1
                                                
        elif j == i+1:
            int_yn[i][j][2][0] = 1
            int_yn[j][i][0][2] = 1
            m = 2
            n = 0
            bond_graph[res2atom_map[i][m]][res2atom_map[j][n]] = 1
            bond_graph[res2atom_map[j][n]][res2atom_map[i][m]] = 1
            
# """Adding disulphide bond"""
# i = 2
# j = 18
# int_yn[i][j][5][5] = 1
# int_yn[j][i][5][5] = 1
# bond_graph[res2atom_map[i][5]][res2atom_map[j][5]] = 1
# bond_graph[res2atom_map[j][5]][res2atom_map[i][5]] = 1

# i = 5
# j = 23
# int_yn[i][j][5][5] = 1
# int_yn[j][i][5][5] = 1
# bond_graph[res2atom_map[i][5]][res2atom_map[j][5]] = 1
# bond_graph[res2atom_map[j][5]][res2atom_map[i][5]] = 1

# i = 9
# j = 25
# int_yn[i][j][5][5] = 1
# int_yn[j][i][5][5] = 1
# bond_graph[res2atom_map[i][5]][res2atom_map[j][5]] = 1
# bond_graph[res2atom_map[j][5]][res2atom_map[i][5]] = 1
            
#%%
"""bond angles"""
count_bonds = 0
for i in range(tot_atoms):
    for j in range(i+1,tot_atoms):
        if bond_graph[i][j] == 1:
            count_bonds += 1

print(count_bonds)  

count_bond_angles = 0

for i in range(tot_atoms):
    for j1 in range(tot_atoms):
        if i == j1 or bond_graph[i][j1] == 0:
            continue
        
        for j2 in range(tot_atoms):
            if i == j2 or j1 == j2 or bond_graph[i][j2] == 0 or bond_graph[j1][j2] == 1:
                continue
            
            res_num1 = atom2res_map[j1][0]
            atom_num1 = atom2res_map[j1][1]
            
            res_num2 = atom2res_map[j2][0]
            atom_num2 = atom2res_map[j2][1]
            
            if int_yn[res_num1][res_num2][atom_num1][atom_num2] == 0 and int_yn[res_num2][res_num1][atom_num2][atom_num1] == 0:
                int_yn[res_num1][res_num2][atom_num1][atom_num2] = 2
                int_yn[res_num2][res_num1][atom_num2][atom_num1] = 2            
                count_bond_angles += 1
                # print(res_num1,atom_num1,res_num2,atom_num2)
        
print(count_bond_angles)   

#%%
"""dihedrals """
count_dihed_angles = 0 

for i in range(tot_atoms):
    for j in range(i,tot_atoms):
        
        if bond_graph[i][j] == 1:
            
            for k1 in range(tot_atoms):
                if j == k1 or bond_graph[i][k1] == 0:
                    continue
   
                for k2 in range(tot_atoms):
                    if i == k2 or bond_graph[j][k2] == 0 or bond_graph[k1][k2] == 1:
                        continue
      
                    res_num1 = atom2res_map[k1][0]
                    atom_num1 = atom2res_map[k1][1]
                    
                    res_num2 = atom2res_map[k2][0]
                    atom_num2 = atom2res_map[k2][1]
                     
                    if int_yn[res_num1][res_num2][atom_num1][atom_num2] == 3 and int_yn[res_num2][res_num1][atom_num2][atom_num1] == 3:
                        count_dihed_angles += 1
                        
                    if int_yn[res_num1][res_num2][atom_num1][atom_num2] == 0 and int_yn[res_num2][res_num1][atom_num2][atom_num1] == 0:
                        int_yn[res_num1][res_num2][atom_num1][atom_num2] = 3
                        int_yn[res_num2][res_num1][atom_num2][atom_num1] = 3
                        count_dihed_angles += 1
                                

print(count_dihed_angles)

                  
            

#%%
allatoms_type = ["" for x in range(tot_atoms)]

vdw_pro = np.zeros((pro_length,pro_length,26,26,3),dtype = float)

for i in range(pro_length):
    aa1 = given_seq[i]
    
    if i == 0:
        aa1 += "NT"
    if i == pro_length - 1:
        aa1 += "CT"
    
    for atom_type1 in list(res_info_all[aa1].keys()):
        m = list(res_info_all[aa1].keys()).index(atom_type1)
        
        atom_num = res2atom_map[i][m]
        allatoms_type[atom_num] = str(res_info_all[aa1][atom_type1]["charmmid"])
        
        
        for j in range(i+1,pro_length):
            aa2 = given_seq[j]
            
            if j == 0:
                aa2 += "NT"
            if j == pro_length - 1:
                aa2 += "CT"
            
            for atom_type2 in list(res_info_all[aa2].keys()):
                n = list(res_info_all[aa2].keys()).index(atom_type2)  
                
                if int_yn[i][j][m][n] in [-1,1,2]:
                    continue
                
                charmmatom1 = str(res_info_all[aa1][atom_type1]["charmmid"])
                charmmatom2 = str(res_info_all[aa2][atom_type2]["charmmid"])
                
                eps1 = vdw_dict[charmmatom1]["eps"]
                sigma1 = vdw_dict[charmmatom1]["rmin"]
                eps2 = vdw_dict[charmmatom2]["eps"]
                sigma2 = vdw_dict[charmmatom2]["rmin"]
                    
                if int_yn[i][j][m][n] == 3:
                    if len(vdw_dict[charmmatom1]) == 4:
                        eps1 = vdw_dict[charmmatom1]["eps14"]
                        sigma1 = vdw_dict[charmmatom1]["rmin14"]
                    if len(vdw_dict[charmmatom2]) == 4:
                        eps2 = vdw_dict[charmmatom2]["eps14"]
                        sigma2 = vdw_dict[charmmatom2]["rmin14"]                    
                    
                vdw_pro[i][j][m][n][0] = np.sqrt(eps1*eps2)
                vdw_pro[i][j][m][n][1] = (sigma1+sigma2)
                vdw_pro[i][j][m][n][2] = pq[i,m]*pq[j,n]                 
    
        
        for atom_type2 in list(res_info_all[aa1].keys()):
            n = list(res_info_all[aa1].keys()).index(atom_type2)  
            if m >= n:
                continue
            
            if int_yn[i][i][m][n] in [-1,1,2]:
                continue    
            
            charmmatom1 = str(res_info_all[aa1][atom_type1]["charmmid"])
            charmmatom2 = str(res_info_all[aa1][atom_type2]["charmmid"])
            
            eps1 = vdw_dict[charmmatom1]["eps"]
            sigma1 = vdw_dict[charmmatom1]["rmin"]
            eps2 = vdw_dict[charmmatom2]["eps"]
            sigma2 = vdw_dict[charmmatom2]["rmin"]
                
            if int_yn[i][i][m][n] == 3:
                if len(vdw_dict[charmmatom1]) == 4:
                    eps1 = vdw_dict[charmmatom1]["eps14"]
                    sigma1 = vdw_dict[charmmatom1]["rmin14"]
                if len(vdw_dict[charmmatom2]) == 4:
                    eps2 = vdw_dict[charmmatom2]["eps14"]
                    sigma2 = vdw_dict[charmmatom2]["rmin14"]                    
                
            vdw_pro[i][i][m][n][0] = np.sqrt(eps1*eps2)
            vdw_pro[i][i][m][n][1] = (sigma1+sigma2)
            vdw_pro[i][i][m][n][2] = pq[i,m]*pq[i,n]             


#%%
tot_energy_CP = 0
tot_energy_VDW = 0
energy_pairs = np.zeros((pro_length,pro_length),dtype = float)

for i in range(pro_length):
    for j in range(i+1,pro_length):
        for m in range(atom_num_pro[i]):
            for n in range(atom_num_pro[j]):
                if int_yn[i][j][m][n] in [-1,1,2,4]:
                    continue
                else:
                    
                    # q1 = pq[i][m]
                    # q2 = pq[j][n]
                    dis = np.linalg.norm(atom_coor[i,m,:] - atom_coor[j,n,:])
                    
                    Energy = 331.21*vdw_pro[i][j][m][n][2]/dis
                    # Energy = 331.21*pq[i][m]*pq[j][n]/dis
                    
                    tot_energy_CP += Energy
                    energy_pairs[i][j] += Energy
                    
                    eps = vdw_pro[i][j][m][n][0]
                    rmin = vdw_pro[i][j][m][n][1]
                    
                    x = (rmin/dis)**6
                    
                    Energy = eps*((x**2) - 2*x)
                    tot_energy_VDW += Energy
                    
                    # if (dis/rmin < 0.80/1.12):
                    #     print(i,j,m,n,dis*1.12/rmin,dis,rmin,Energy)
                    
                    energy_pairs[i][j] += Energy
        
    
    for m in range(atom_num_pro[i]):
        for n in range(m+1,atom_num_pro[i]):
            if int_yn[i][i][m][n] in [-1,1,2,4]:
                continue
            else:
                # q1 = pq[i][m]
                # q2 = pq[i][n]
                dis = np.linalg.norm(atom_coor[i,m,:] - atom_coor[i,n,:])
                
                Energy = 331.21*vdw_pro[i][i][m][n][2]/dis
                # Energy = 331.21*pq[i][m]*pq[i][n]/dis
                  
                tot_energy_CP += Energy
                energy_pairs[i][i] += Energy
                
                eps = vdw_pro[i][i][m][n][0]
                rmin = vdw_pro[i][i][m][n][1]
                
                x = (rmin/dis)**6
                
                Energy = eps*((x**2) - 2*x)
                energy_pairs[i][i] += Energy
                tot_energy_VDW += Energy                
 
#%%
os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = seq_name+".a2r"

with open(filename, 'w') as fh:
    fh.write('{:d}\n'.format(tot_atoms))
    for i in range(tot_atoms):
        fh.write('{:d} {:d} {:d}\n'.format(i,atom2res_map[i][0],atom2res_map[i][1]))
        
#%%       
os.chdir(os.path.dirname(full_path))
os.chdir('../files4RRR')

filename = seq_name+".par" 

with open(filename, 'w') as fh:
    for i in range(pro_length):
        aa1 = given_seq[i]        
        if i == 0:
            aa1 += "NT"
        if i == pro_length - 1:
            aa1 += "CT"
            
        for atom_type1 in list(res_info_all[aa1].keys()):
            
            m = list(res_info_all[aa1].keys()).index(atom_type1)

            charmmatom1 = str(res_info_all[aa1][atom_type1]["charmmid"])          
            eps = vdw_dict[charmmatom1]["eps"]
            sigma = vdw_dict[charmmatom1]["rmin"]            
            
            fh.write('{:d} {:d} {:.6f} {:.6f} {:.6f}\n'.format(i,m,pq[i][m],eps,sigma))

#%%       
os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = seq_name+".nonb" 

with open(filename, 'w') as fh:
    for i in range(pro_length):
        for m in range(atom_num_pro[i]):
            for n in range(m+1,atom_num_pro[i]):
                if int_yn[i][i][m][n] in [0,3]:
                    fh.write('{:d} {:d} {:.6f} {:.6f} {:.6f}\n'.format(res2atom_map[i][m],res2atom_map[i][n],vdw_pro[i][i][m][n][2],vdw_pro[i][i][m][n][0],vdw_pro[i][i][m][n][1]))
                
        
        for j in range(i+1,pro_length):
            for m in range(atom_num_pro[i]):
                for n in range(atom_num_pro[j]):
                    if int_yn[i][j][m][n] in [0,3]:
                        fh.write('{:d} {:d} {:.6f} {:.6f} {:.6f}\n'.format(res2atom_map[i][m],res2atom_map[j][n],vdw_pro[i][j][m][n][2],vdw_pro[i][j][m][n][0],vdw_pro[i][j][m][n][1]))
                        
  
#%%
"""creating the CT and SC Motif library"""

CT_motifs = np.zeros((20,50,6,3),dtype = float)
CT_count = np.zeros(20,dtype = int)

max_motifs = 20
SC_motifs = np.zeros((20,np.max(numF),max_motifs,np.max(sizeF),3),dtype = float)
count_SC_motifs = np.zeros(20,dtype = int)

protein_list = ["2JOF","1LVZ","4FC1","1LB0"]

os.chdir(os.path.dirname(full_path))
os.chdir('../PDB')

for protein in protein_list:
    filename = protein+".pdb"  

    mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
    parser = PDBParser()
    structure = parser.get_structure(protein, filename)  
    model = structure[0]
    # dssp = DSSP(model, filename, dssp='mkdssp')    # for the dssp module
    
    id_chain = list(model.child_dict.keys())[0]
    chain = model[id_chain]


    all_atoms = Selection.unfold_entities(structure, "A")
    all_residues = Selection.unfold_entities(chain, "R")
    
    ppb=PPBuilder()
    
    seq = ""
    for pp in ppb.build_peptides(structure):
        seq += str(pp.get_sequence())
    
    num_res = 0
    
    for residue in all_residues:
        if (residue.get_parent()).get_id() == id_chain and residue.get_id()[0] == ' ':
            num_res += 1


    temp_coor = np.zeros((num_res,26,3),dtype = float)
    temp_num_atoms = np.zeros(num_res,dtype = int)

    for i in range(num_res):
        res = all_residues[i]
        if (res.get_parent()).get_id() != id_chain:
            continue
    
        if res.get_id()[0] != ' ':
            continue  
        
        a_key = (id_chain,res.get_id())
        aa = seq[i]
        
        temp_num_atoms[i] = atom_num_AA[residue_type(aa)]
    
        if i == 0:
            aa += "NT"
            temp_num_atoms[i] += 2
        if i == num_res - 1:
            aa += "CT"
            temp_num_atoms[i] += 1
    
        for atom_type in list(res_info_all[aa].keys()):
            j = list(res_info_all[aa].keys()).index(atom_type)
            
            if atom_type[0] == "H":
                try:
                    temp_coor[i,j,:] = res[atom_type].get_coord()
                except:
                    atom_type_new = atom_type.replace("H", "D", 1)
                    temp_coor[i,j,:] = res[atom_type_new].get_coord()
                    
            else:
                temp_coor[i,j,:] = res[atom_type].get_coord()


    for i in range(num_res-1):
        aa = residue_type(seq[i+1])
        
        if aa == 10:
            sizeCT = 5
        else:
            sizeCT = 6
    
        coor = np.zeros((sizeCT,3),dtype = float)
        
        coor[0,:] = temp_coor[i,1,:]
        coor[1,:] = temp_coor[i,2,:]
        coor[2,:] = temp_coor[i,3,:]
        coor[3,:] = temp_coor[i+1,0,:]
        coor[4,:] = temp_coor[i+1,1,:]
        
        if sizeCT == 6:
            coor[5,:] = temp_coor[i+1,atom_num_AA[aa]-1,:]
    
        coor = trans_origin(coor)
        
        CT_motifs[aa][CT_count[aa]][0:sizeCT,:] = coor[:,:] 
        CT_count[aa] += 1


    for i in range(num_res):
        aa = seq[i]
        a_SC = residue_type(aa)
    
        try:
            for frags in range(0,numF[a_SC]):
                num_atoms = sizeF[a_SC][frags]
                coord = np.zeros((num_atoms,3),dtype = float)
                
                atom_no = 0
                motif_atom = 0
                
                for atom_no in atom_numF[a_SC][frags][0:sizeF[a_SC][frags]]:
                    coord[motif_atom,:] = temp_coor[i,atom_no,:]
                    motif_atom += 1
                    
                if motif_atom != num_atoms:
                    continue
                
                coord[:,:] = trans_origin(coord)
                
                SC_motifs[a_SC,frags,count_SC_motifs[a_SC],0:num_atoms,:] = coord[:,:] 
    
            count_SC_motifs[a_SC] += 1
          
        except:
            continue 
#%%
os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = "CTMotifs.txt"

with open(filename, 'w') as fh:
    for i in range(20):        
        
        if i == 10:
            sizeCT = 5
        else:
            sizeCT = 6
                        
        if CT_count[i] == 0:
            continue

        fh.write('{} {:d}\n'.format(AA_list[i],CT_count[i])) 
        for z in range(CT_count[i]):

            for q in range(sizeCT):
                fh.write('{:.6f} \t {:.6f} \t {:.6f}\n'.format(CT_motifs[i,z,q,0],CT_motifs[i,z,q,1],CT_motifs[i,z,q,2]))                   
            fh.write('\n')  
            


#%%

"""writing SC Motif library"""

os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

AA_list = "RNDQEHKCMGPSTWYILVFA"

filename = "SCMotifs.txt"

with open(filename, 'w') as fh:
    for i in range(20):
        for frags in range(0,numF[i]):
            num_atoms = sizeF[i][frags]         
                        
            if count_SC_motifs[i] == 0:
                continue
            
            # fh.write('{} {:d} {}\n'.format(AA_list[i],frags,len(np.where(motifs_state[i][frags] == 2)[0])))
            fh.write('{} {:d} {}\n'.format(AA_list[i],frags,count_SC_motifs[i])) 
            for z in range(count_SC_motifs[i]):
                # if motifs_state[i][frags][z] == 2:
                                        
                for q in range(num_atoms):
                    fh.write('{:.6f} \t {:.6f} \t {:.6f}\n'.format(SC_motifs[i,frags,z,q,0],SC_motifs[i,frags,z,q,1],SC_motifs[i,frags,z,q,2]))                   
                fh.write('\n')  
                

#%%

"""Dihedral information"""

os.chdir(os.path.dirname(full_path))
os.chdir('../GenFiles')


from parmed.charmm import CharmmParameterSet
# from parmed.charmm import read_topology_file
# from parmed.charmm import load_set
# params = CharmmParameterSet('par_all22_prot.inp', 'top_all22_prot.inp')
# params = CharmmParameterSet.load_set()
prms = CharmmParameterSet('par_all36m_prot.prm')  
 

#%%
"""dihedrals """
# count_dihed_angles = 1000
dihed_coors = np.zeros((count_dihed_angles,4,3),dtype = float)
dihed_atoms = np.zeros((count_dihed_angles,4),dtype = int)
dihed_prms = np.zeros((count_dihed_angles,5,3),dtype = float)
total_dihed_prms = np.zeros(count_dihed_angles,dtype = int)

count_dh = 0
for i in range(tot_atoms):
    for j in range(i,tot_atoms):
        
        if bond_graph[i][j] == 1:
            
            for k1 in range(tot_atoms):
                if j == k1 or bond_graph[i][k1] == 0:
                    continue
   
                for k2 in range(tot_atoms):
                    if i == k2 or bond_graph[j][k2] == 0 or bond_graph[k1][k2] == 1:
                        continue
      
                    res_num1 = atom2res_map[k1][0]
                    atom_num1 = atom2res_map[k1][1]
                    
                    res_num2 = atom2res_map[k2][0]
                    atom_num2 = atom2res_map[k2][1]
                    
                    if int_yn[res_num1][res_num2][atom_num1][atom_num2] == 3 and int_yn[res_num2][res_num1][atom_num2][atom_num1] == 3:                        
                        dihed_atoms[count_dh][0] = k1 
                        dihed_atoms[count_dh][1] = i
                        dihed_atoms[count_dh][2] = j
                        dihed_atoms[count_dh][3] = k2
                        
                        at1 = allatoms_type[k1]
                        at2 = allatoms_type[i]
                        at3 = allatoms_type[j]
                        at4 = allatoms_type[k2]
                        
                        try:
                            total_dihed_prms[count_dh] = len(prms.dihedral_types[(at1, at2, at3, at4)])
                            
                            for m in range(total_dihed_prms[count_dh]):
                                dihed_prms[count_dh][m][0] = prms.dihedral_types[(at1, at2, at3, at4)][m].phi_k   
                                dihed_prms[count_dh][m][1] = prms.dihedral_types[(at1, at2, at3, at4)][m].per
                                dihed_prms[count_dh][m][2] = prms.dihedral_types[(at1, at2, at3, at4)][m].phase
                            
                            if total_dihed_prms[count_dh] > 1:
                                print(at1,at2,at3,at4)
                                print(prms.dihedral_types[(at1, at2, at3, at4)])
                                
                        except:
                            # print(at1,at2,at3,at4, k1,i,j,k2)
                            # print(at1,at2,at3,at4,"X")
                            at1 = 'X'
                            at4 = 'X'
                            
                        try:
                            if at1 == 'X' and at4 == 'X':
                                total_dihed_prms[count_dh] = len(prms.dihedral_types[(at1, at2, at3, at4)])
                                print(total_dihed_prms[count_dh])
                                for m in range(total_dihed_prms[count_dh]):
                                    dihed_prms[count_dh][m][0] = prms.dihedral_types[(at1, at2, at3, at4)][m].phi_k   
                                    dihed_prms[count_dh][m][1] = prms.dihedral_types[(at1, at2, at3, at4)][m].per
                                    dihed_prms[count_dh][m][2] = prms.dihedral_types[(at1, at2, at3, at4)][m].phase

                        except:
                            print(allatoms_type[k1],at2,at3,allatoms_type[k2])                            
                   
                        if dihed_prms[count_dh][0][0] < 0:
                            print(total_dihed_prms[count_dh],at1,at2,at3,at4)
                        count_dh += 1
         
#%%
def dihed_calc(vec1,vec2,vec3,vec4):
    vecBA = vec1 - vec2
    vecCB = vec2 - vec3
    vecCD = vec4 - vec3
    
    perp_ABC = np.cross(vecBA, vecCB)
    perp_BCD = np.cross(vecCD, vecCB)
    
    dot = np.dot(perp_ABC,perp_BCD)
    dt = -np.linalg.det(np.dstack([vecCB,perp_ABC,perp_BCD]))
    dt = dt/np.linalg.norm(vecCB)
    
    return math.atan2(dt,dot)

#%%
tot_energy_dihed = 0
tkb = 0
for i in range(count_dihed_angles):
    v1 = dihed_coors[i][0]
    v2 = dihed_coors[i][1]
    v3 = dihed_coors[i][2]
    v4 = dihed_coors[i][3]
    
    chim = dihed_calc(v1,v2,v3,v4)
    
    v1 = Vector(dihed_coors[i][0])
    v2 = Vector(dihed_coors[i][1])
    v3 = Vector(dihed_coors[i][2])
    v4 = Vector(dihed_coors[i][3])
    
    at1 = allatoms_type[dihed_atoms[i][0]]
    at2 = allatoms_type[dihed_atoms[i][1]]
    at3 = allatoms_type[dihed_atoms[i][2]]
    at4 = allatoms_type[dihed_atoms[i][3]]
    
    chi = calc_dihedral(v1, v2, v3, v4) 
    
    print(i,chi,chim)
    
    kchi = dihed_prms[i][0]
    n = dihed_prms[i][1]
    delta = dihed_prms[i][2]
    
    energy_dihed = kchi*(1+np.cos((n*chi) - math.radians(delta)))
    tot_energy_dihed += energy_dihed
    
    tkb += kchi
    
    # print("{:d} {} {} {} {} {:0.2f} {:0.2f} {:0.2f}\n".format(i,at1,at2,at3,at4,chi,energy_dihed,kchi))


 
#%%
"""checking bond energies"""  

bond_atoms = np.zeros((count_bonds,2),dtype = int)
bond_prms = np.zeros((count_bonds,2),dtype = float)
tot_energy_bond = 0
z = 0

for i in range(tot_atoms):
    for j in range(i+1,tot_atoms):
        if bond_graph[i][j] == 1:
            
            bond_atoms[z][0] = i
            bond_atoms[z][1] = j
            
            res_num1 = atom2res_map[i][0]
            atom_num1 = atom2res_map[i][1]
            
            res_num2 = atom2res_map[j][0]
            atom_num2 = atom2res_map[j][1]
            
            at2 = allatoms_type[i]
            at3 = allatoms_type[j]
            
            # bond_length = np.linalg.norm(atom_coor[res_num1][atom_num1] - atom_coor[res_num2][atom_num2])
            
            kb = prms.bond_types[(at2, at3)].k
            bo = prms.bond_types[(at2, at3)].req
            
            bond_prms[z][0] = kb
            bond_prms[z][1] = bo            
            
            # energy_bond = kb*((bond_length - bo)**2)
            # tot_energy_bond += energy_bond
            
            z += 1
            
            
            
            # print("{:d} {} {} {:0.2f} {:0.2f} {:0.2f} {:0.2f}\n".format(i,at2,at3,kb,bo,bond_length,energy_bond))

#%%    
os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = seq_name+".dhd"

with open(filename, 'w') as fh:
    
    fh.write('{:d}\n'.format(count_bonds)) 
    
    for i in range(count_bonds):
        fh.write('{:d} {:d} {:.6f} {:.6f} \n'.format(bond_atoms[i][0],bond_atoms[i][1],bond_prms[i][0],bond_prms[i][1]))
    
    
    fh.write('{:d}\n'.format(count_dihed_angles)) 
    
    for i in range(count_dihed_angles): 
        fh.write('{:d} \t {:d} {:d} {:d} {:d} \n'.format(total_dihed_prms[i],dihed_atoms[i][0],dihed_atoms[i][1],dihed_atoms[i][2],dihed_atoms[i][3]))
        for m in range(total_dihed_prms[i]):
            fh.write('{:.6f} {:.6f} {:.6f}\n'.format(dihed_prms[i][m][0],dihed_prms[i][m][1],math.radians(dihed_prms[i][m][2])))
 

#%%
"""Coordinate File"""

os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = seq_name+".coor" 

with open(filename, 'w') as fh:
    for i in range(pro_length):
        for m in range(atom_num_pro[i]):
            fh.write('{:.6f} {:.6f} {:.6f}\n'.format(atom_coor[i,m,0],atom_coor[i,m,1],atom_coor[i,m,2]))
        fh.write('\n')

    # for i in range(num_sol):
    #     for j in range(3):
    #         fh.write('{:.6f} {:.6f} {:.6f}\n'.format(sol_coor[i,j,0],sol_coor[i,j,1],sol_coor[i,j,2]))

        # fh.write('\n') 
#%%
"""Sequence File"""  
num_sol = 0 
os.chdir(os.path.dirname(full_path))
os.chdir('../Files4RRR')

filename = seq_name+".seq"

with open(filename, 'w') as fh:
    fh.write('{:d} {:d}\n'.format(pro_length,num_sol)) 
    for i in range(pro_length):
        fh.write('{}'.format(given_seq[i])) 
    fh.write('\n') 