import re,os,sys,glob
import numpy as np
from scipy.stats import mannwhitneyu,pearsonr,ttest_ind
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns 

baseComplement = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }
def revc(seq):
    return "".join([baseComplement[base] for base in seq[::-1]])

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]   
        return Data


def IR_function(sequence,min_spacer_size,max_spacer_size,arm_size):
	# Function to find IRs
        found=False
        positions=[]    
        for spacer in range(min_spacer_size,max_spacer_size,1):
                size = spacer+arm_size*2
                for k in range(len(sequence)-size):
                        seq_part1=sequence[k:k+arm_size]
                        seq_part2=revc(sequence[k+arm_size+spacer:k+spacer+2*arm_size])
                        if seq_part1==seq_part2:
                                        found=True
                                        positions+=[k]
        return found#,positions


def DR_function(sequence,min_spacer_size,max_spacer_size,arm_size):
	# Function to find DRs
        found=False
        positions=[]    
        for spacer in range(min_spacer_size,max_spacer_size,1):
                size = spacer+arm_size*2
                for k in range(len(sequence)-size):
                        seq_part1=sequence[k:k+arm_size]
                        seq_part2=sequence[k+arm_size+spacer:k+spacer+2*arm_size]
                        if seq_part1==seq_part2 and find_reps(seq_part1)<0.8:
                                #print(seq_part1,seq_part2,
                                found=True
                                positions+=[k]
        return found#,positions

def function_MR(sequence,min_spacer_size,max_spacer_size,arm_size):
        found=False
        positions=[]
        for spacer in range(min_spacer_size,max_spacer_size,1):
                size = spacer+arm_size*2
                for k in range(len(sequence)-size):
                        seq_part1=sequence[k:k+arm_size]
                        seq_part2=sequence[k+arm_size+spacer:k+spacer+2*arm_size][::-1]
                        if seq_part1==seq_part2:
                                found=True
                                positions+=[k]
        return found#,positions

def function_STR(sequence):
        WithSTR=False
        for size in range(1,10):
                distance=0
                pos = 0
                motif = sequence[pos:pos+size]
                repeat = sequence[pos:pos+size]
                counter=0
                while distance<200:
                        while repeat==motif:
                                counter+=1
                                pos=pos+len(repeat)
                                repeat = sequence[pos:pos+size]
                                distance+=len(repeat)
                        else:
                                distance+=len(repeat)
                                if counter>=5 and len(motif)*counter>=12:
                                        #print repeat, counter
                                        WithSTR=True
                                counter=0
                                motif = sequence[pos:pos+size]
                                repeat = sequence[pos:pos+size]
        return WithSTR


def find_reps(seq):
        countsL=[]
        for s in range(1,5):
                repeat =seq[:s]
                counts=0
                for v in range(0,len(seq),len(repeat)):
                        if seq[v:v+len(repeat)]==repeat:
                                counts+=1
                countsL+=[(s*counts)/float(len(seq))]

        return max(countsL)

def Z_DNA_find(seq):
        scores=[]
        purines=["A","G"]
        pyrimidines=["T","C"]
        score=1
        for k in range(1,len(seq),1):
                if ((seq[k-1] in purines and seq[k] in pyrimidines) or (seq[k-1] in pyrimidines and seq[k] in purines)) and seq[k-1:k+1]!="AT" and seq[k-1:k+1]!="TA":
                        score+=1
                        if score==len(seq)-1:
                                scores+=[score]
                else:
                        scores+=[score]
                        score=1
        if max(scores)>=10:
                return True
        else:
                return False
