import re,os,sys
import glob
from random import shuffle
import random,pybedtools


def sequence_from_coordinates(coordinates):
        path_file="hg19/"
        chromosome_to_use=coordinates[0]
        chromosome_seq_=read_file(path_file+"/"+str(chromosome_to_use)+".fa")
        sequence=chromosome_seq_[coordinates[1]:coordinates[2]]
        return sequence

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def read_file(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Seqs=[]
        for i in data[1:]:
                Seqs+=[i.strip()]
        Seq="".join(Seqs)
        return Seq.upper()

def in_region(lineR,DataRegion):
        found=False
        for k in DataRegion:
                if k[0]==lineR[0]:
                        if int(lineR[1])>=int(k[1]) and int(lineR[1])<int(k[2]):
                                if int(lineR[2])>=int(k[1]) and int(lineR[2])<int(k[2]):
                                        found=True
        return found

Data_include=reader("hg19-v0-wgs_evaluation_regions.v1.interval_list.bed")

files=glob.glob("gnomad.genomes.r2.1.1.sites_insertions.bed.PASS")
for one in file:
        DataL=reader(one)
        parts=30000
        total=len(DataL)
        part=(int(sys.argv[1])-1)
        step=int(total*float(1/float(parts)))
        RepeatsL=DataL[part*step:(part+1)*step]
        if part==parts-1:
                RepeatsL=DataL[part*step:]
        if part>parts or part<0:
                break

        datafile=open("controls_gnomads/"+one.split("/")[-1].split("_")[1]+"_control_"+str(part),"w")
        DataL=RepeatsL
        for lined in DataL:
                seq = sequence_from_coordinates([lined[0],int(lined[1])-100,int(lined[2])+100]).upper()
                GC_content = (seq.count("G")+seq.count("C"))/float(len(seq))
                size=int(lined[2])-int(lined[1])
                sequence_F=None
                seq_coordinates=None
                score=False
                GC_content_C=-2
                window=10000
                while abs(GC_content_C-GC_content)<=0.025 or score!=True:
                        direction =random.choice([0,1])
                        if direction==0:
                                coordinates = range(int(lined[1])-window,int(lined[1])-3)
                                chosen_pos  = random.choice(coordinates)
                                seq_coordinates = [lined[0],chosen_pos,chosen_pos+size]
                                seq_control=sequence_from_coordinates([seq_coordinates[0],int(seq_coordinates[1])-100,int(seq_coordinates[2])+100]).upper()
                                GC_content_C = (seq_control.count("G")+seq_control.count("C"))/float(len(seq_control))
                                score=in_region(seq_coordinates,Data_include)

                        elif direction==1:
                                coordinates = range(int(lined[1])+3,int(lined[1])+window)
                                chosen_pos  = random.choice(coordinates)
                                seq_coordinates = [lined[0],chosen_pos,chosen_pos+size]
                                seq_control=sequence_from_coordinates([seq_coordinates[0],int(seq_coordinates[1])-100,int(seq_coordinates[2])+100]).upper()
                                GC_content_C = (seq_control.count("G")+seq_control.count("C"))/float(len(seq_control))
                                score=in_region(seq_coordinates,Data_include)
                else:
                        datafile.write(str(seq_coordinates[0])+'\t'+str(int(seq_coordinates[1]))+'\t'+str(int(seq_coordinates[2]))+'\n')
        datafile.close()

