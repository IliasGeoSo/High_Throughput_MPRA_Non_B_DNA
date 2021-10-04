import re,os,sys
import glob
from random import shuffle
import random,pybedtools

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
                                found=True
        return found

Data_include=reader("hg19-v0-wgs_evaluation_regions.v1.interval_list.bed")
part=(int(sys.argv[1])-1)
datafile=open("gnomad_control_"+str(part),"w")
files=glob.glob("chroms/*")
for one in files:
        DataL=reader(one)
        parts=45000
        total=len(DataL)
        step=int(total*float(1/float(parts)))
        RepeatsL=DataL[part*step:(part+1)*step]
        if part==parts-1:
                RepeatsL=DataL[part*step:]
        if part>parts or part<0:
                break

        DataL=RepeatsL

        path_file="hg19/"
        chromosome_to_use=DataL[0][0]
        chromosome_seq_=read_file(path_file+"/"+str(chromosome_to_use)+".fa")

        for lined in DataL:
                seq = chromosome_seq_[int(lined[1])-1:int(lined[2])+1]
                sequence_F=None
                seq_coordinates=None
                score=False
                window=10000
                while seq != sequence_F or score!=True:
                        direction =random.choice([0,1])
                        if direction==0:
                                coordinates = range(int(lined[1])-window,int(lined[1])-3)
                                chosen_pos  = random.choice(coordinates)
                                seq_coordinates = [lined[0],chosen_pos-1,chosen_pos+2]
                                sequence_F=chromosome_seq_[chosen_pos-1:chosen_pos+2]
                                score=in_region(seq_coordinates,Data_include)

                        elif direction==1:
                                coordinates = range(int(lined[1])+3,int(lined[1])+window)
                                chosen_pos  = random.choice(coordinates)
                                seq_coordinates = [lined[0],chosen_pos-1,chosen_pos+2]
                                sequence_F=chromosome_seq_[chosen_pos-1:chosen_pos+2]
                                score=in_region(seq_coordinates,Data_include)
                else:
                        datafile.write(str(seq_coordinates[0])+'\t'+str(int(seq_coordinates[1])+1)+'\t'+str(int(seq_coordinates[2])-1)+'\t'+seq+'\n')
datafile.close()

