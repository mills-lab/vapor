from __future__ import print_function

import os
import sys
import itertools,math,random,numpy,scipy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
#import vapor_vali.plotting as plotting
#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
from scipy.cluster.vq import vq, kmeans, whiten
from scipy import stats
from scipy.stats import linregress
from scipy.spatial import distance
from sklearn import cluster
global invert_base
invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C','N' : 'N','a' : 't', 't' : 'a', 'c' : 'g', 'g' : 'c','n' : 'n'}
global default_flank_length
default_flank_length=500
global default_read_length
default_read_length=4000    #average length of pacbio read
global default_max_sv_test
default_max_sv_test=10000 #when size of a sv block excessed default_max_sv_test, try junctions instead of event

def alt_seq_readin(ref,info,flank_length):
    #eg of info=('abcde/abcde', 'db/db', '1', '14436355', '14436655', '14437075', '14438241', '14438566', '14438940']
    block_hash={}
    for x in range(len(info)-4):
        block_hash[chr(97+x)]=[info[2],int(info[x+3]),int(info[x+4])]
    block_hash['-']=[block_hash[info[0][0]][0],block_hash[info[0][0]][1]-flank_length,block_hash[info[0][0]][1]]
    block_hash['+']=[block_hash[info[0][-1]][0],block_hash[info[0][-1]][2],block_hash[info[0][-1]][2]+flank_length]
    if not info[1].split('/')[0]==info[1].split('/')[1]:
        out=[ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2]),ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2])]
        rec=-1
        for x in info[1].split('/'):
            rec+=1
            temp_x=[]
            for y in x:
                if not y=='^':temp_x.append(y)
                else:    temp_x[-1]+=y
            for y in temp_x:
                if not '^' in y:
                    out[rec]+=ref_seq_readin(ref,block_hash[y][0],block_hash[y][1],block_hash[y][2])
                else:
                    out[rec]+=ref_seq_readin(ref,block_hash[y[0]][0],block_hash[y[0]][1],block_hash[y[0]][2],'TRUE')
        out[0]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        out[1]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        return out
    else:
        out=[ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2])]
        rec=-1
        for x in [info[1].split('/')[0]]:
            rec+=1
            temp_x=[]
            for y in x:
                if not y=='^':temp_x.append(y)
                else:    temp_x[-1]+=y
            for y in temp_x:
                if not '^' in y:
                    out[rec]+=ref_seq_readin(ref,block_hash[y][0],block_hash[y][1],block_hash[y][2])
                else:
                    out[rec]+=ref_seq_readin(ref,block_hash[y[0]][0],block_hash[y[0]][1],block_hash[y[0]][2],'TRUE')
        out[0]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        return [out[0],out[0]]

def bam_in_decide(bam_in,bps):
    if os.path.isfile(bam_in):
        return [bam_in]
    else:
        temp_bam_in=[]
        bam_in_path='/'.join(bam_in.split('/')[:-1])+'/'
        if 'XXX' in bam_in.split('/')[-1]:
            bam_in_keys=bam_in.split('/')[-1].split('XXX')
        elif '*' in bam_in.split('/')[-1]:
            bam_in_keys=bam_in.split('/')[-1].split('*')
        else:
            print ('Error: invalid name for pacbio files !')
        for k1 in os.listdir(bam_in_path):
            if k1.split('.')[-1]==bam_in.split('.')[-1]:
                flag=0
                for y in bam_in_keys:
                    if not y in k1:
                        flag+=1
                if flag==0:
                    temp_bam_in.append(bam_in_path+k1)
        return temp_bam_in

def block_around_check(alt_allele,ref_allele):
    #eg of alt_allele='abcab'  eg of ref_allele='abcd'
    alt_juncs=[[['-']+letter_split(alt_allele)+['+']][0][j:j+2] for j in range(len(letter_split(alt_allele))+1)]
    ref_juncs=[[['-']+letter_split(ref_allele)+['+']][0][j:j+2] for j in range(len(letter_split(alt_allele))+1)]
    new_juncs=[i for i in alt_juncs if not i in ref_juncs]
    return new_juncs

def bp_to_chr_hash(bps,chromos,flank_length=500):
    #eg of bps=['chr16', '34910548', '34911339', '34913149', '34913438', '36181068', '36181482']
    temp1=[]
    for i in bps:
        if i in chromos:
            temp1.append([i])
        else:
            temp1[-1].append(i)
    out={}
    rec=-1
    for k1 in temp1:
        for k2 in range(len(k1[2:])):
            rec+=1
            out[chr(97+rec)]=[k1[0],k1[k2+1],k1[k2+2]]
    out['+']=[out[sorted(out.keys())[-1]][0],out[sorted(out.keys())[-1]][2],str(int(out[sorted(out.keys())[-1]][2])+flank_length)]
    out['-']=[out['a'][0],str(int(out['a'][1])-flank_length),int(out['a'][1])]
    return out

def bp_to_block_len(bp_info):
    #eg of bp_info=['chr1',127906515, 127907283, 127907399]
    out={}
    for i in range(len(bp_info)-2):
        out[chr(97+i)]=bp_info[i+2]-bp_info[i+1]
    return out

def block_modify(block,chromos):
    #eg of block=['chr16', '34911339', '34913149', 'chr16', '34913149', '34913438']
    out=[]
    for x in block:
        if x in chromos:
            if out==[]:    out.append([x])
            else:
                if not x in out[-1]:    out.append([x])
        else:   out[-1].append(x)
    out_new=[]
    for x in out:
        out_new.append([])
        for y in x:
            if x.count(y)==1:
                out_new[-1].append(y)
    out_new_2=[]
    for x in out_new:
        if len(x)==3:
            out_new_2.append(x)
        else:
            for y in range(int((len(x)-1)/2)):
                out_new_2.append([x[0],x[2*y+1],x[2*y+2]])
    return out_new_2

def block_subsplot(bp_list,chromos):
    #eg of bp_list=['chr1', '9061390', '9061470', 'chr14', '93246136', '93248034']
    out=[]
    for x in bp_list:
        if not x in chromos:    out[-1].append(int(x))
        else:   out.append([x])
    return out

def calcu_log10(x):
    if x==0:
        return 0
    else:
        return np.log10(x)

def calcu_vapor_single_read_score_abs_dis_m1(ref_seq,alt_seq,x,window_size):
    ref_seq=ref_seq.upper()
    alt_seq=alt_seq.upper()
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_m1(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_m1(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size):
    ref_seq=ref_seq.upper()
    alt_seq=alt_seq.upper()
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if len(ref_dotdata)>2 and len(alt_dotdata)>2:
        if float(len(ref_dotdata))/min([float(len(ref_seq)),float(len(alt_seq))])>0.1: 
            if float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.6 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.6:
                ref_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(ref_dotdata)
                alt_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(alt_dotdata)
                if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
                    return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
                else:    
                    return [0,0]
            else:
                if float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.6: return [1.1,2.1]
                elif float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.6: return [2.1,1.1]
                else: return [0,0]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_directed_dis_m1b(ref_seq,alt_seq,x,window_size):
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        #[ref_clean_dotdata,ref_kept_segs]=clean_dotdata_diagnal_m1b(ref_dotdata)
        #[alt_clean_dotdata,alt_kept_segs]=clean_dotdata_diagnal_m1b(alt_dotdata)
        #ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        #alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        #[ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        #[alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        #ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        #alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        ref_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(ref_dotdata)
        alt_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_reg_calcu(ref_clean_dotdata),eu_dis_reg_calcu(alt_clean_dotdata)]
            #return [eu_dis_abs_calcu(take_off_symmetric_dots(ref_clean_dotdata)),eu_dis_abs_calcu(take_off_symmetric_dots(alt_clean_dotdata))]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_directed_dis_m1b_not_really(ref_seq,alt_seq,x,window_size,ref_bps,alt_bps):
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        ref_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(ref_dotdata)
        alt_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [abs(eu_dis_region_calcu(ref_clean_dotdata,ref_bps)),abs(eu_dis_region_calcu(alt_clean_dotdata,alt_bps))]
            #return [eu_dis_abs_calcu(take_off_symmetric_dots(ref_clean_dotdata)),eu_dis_abs_calcu(take_off_symmetric_dots(alt_clean_dotdata))]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_directed_dis_m1b_redefine_diagnal(ref_seq,alt_seq,x,window_size):
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        ref_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(ref_dotdata)
        alt_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            ref_new_intercept=dis_to_diagnal_most_abundant_defined(ref_clean_dotdata)
            alt_new_intercept=dis_to_diagnal_most_abundant_defined(alt_clean_dotdata)
            return [abs(eu_dis_dir_calcu([[i[0]+ref_new_intercept,i[1]] for i in ref_clean_dotdata])),
                    abs(eu_dis_dir_calcu([[i[0]+alt_new_intercept,i[1]] for i in alt_clean_dotdata]))]
            #return [abs(eu_dis_region_calcu(ref_clean_dotdata,ref_bps)),abs(eu_dis_region_calcu(alt_clean_dotdata,alt_bps))]
            #return [eu_dis_abs_calcu(take_off_symmetric_dots(ref_clean_dotdata)),eu_dis_abs_calcu(take_off_symmetric_dots(alt_clean_dotdata))]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_directed_dis_m1b_maybe(ref_seq,alt_seq,x,window_size,dup_block_bps):
    #ref_ref_dotdata=dotdata(window_size,ref_seq[x[1]:],ref_seq[x[1]:])
    #alt_alt_dotdata=dotdata(window_size,alt_seq[x[1]:],alt_seq[x[1]:])
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    #ref_ref_lines=ref_ref_deviate_lines_describe(ref_ref_dotdata)
    #alt_alt_lines=ref_ref_deviate_lines_describe(alt_alt_dotdata)
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        ref_clean_dotdata=(ref_dotdata)
        alt_clean_dotdata=clean_dotdata_diagnal_and_anti_diagnal(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [abs(eu_dis_abs_calcu(ref_clean_dotdata)),abs(eu_dis_reg_dup_block_calcu(alt_clean_dotdata,dup_block_bps))]
            #return [eu_dis_abs_calcu(take_off_symmetric_dots(ref_clean_dotdata)),eu_dis_abs_calcu(take_off_symmetric_dots(alt_clean_dotdata))]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size):
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if max([float(len(ref_dotdata))/float(len(ref_seq)),float(len(alt_dotdata))/float(len(alt_seq))])>0.1:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_diagnal_m1b(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_diagnal_m1b(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_dots_within_10perc(alt_clean_dotdata),eu_dis_dots_within_10perc(ref_clean_dotdata)]    #if prediction correct [large, small]
        else:    
            return [0,0]
    else:
        return [0,0]

def calcu_vapor_single_read_score_abs_dis_m2(ref_seq,alt_seq,x,window_size):
    ref_dotdata=dotdata(window_size,x[0],ref_seq[x[1]:])
    alt_dotdata=dotdata(window_size,x[0],alt_seq[x[1]:])
    if float(len(ref_dotdata))/float(len(ref_seq))>0.1 and float(len(alt_dotdata))/float(len(alt_seq))>0.1 and float(ref_dotdata[-1][0]-ref_dotdata[0][0])/float(len(ref_seq))>0.7 and float(alt_dotdata[-1][0]-alt_dotdata[0][0])/float(len(alt_seq))>0.7:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_m2(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_m2(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
        else:    
            return [0,0]
    else:
        return [0,0]

def cigar2alignstart_by_pos(cigar,align_start,start,end):
    #eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
    import re
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    read_rec=0
    align_rec=align_start
    cigar_record=''
    for x in cigars:
        if x[1]=='S':
            read_rec+=int(x[0])
        if x[1] in ['M','=']:
            read_rec+=int(x[0])
            align_rec+=int(x[0])
        if x[1]=='D':
            align_rec+=int(x[0])
        if x[1]=='I':
            read_rec+=int(x[0])
        cigar_record=x
        if align_rec>start-1: break
    start_dis=int(align_rec)-start
    if cigar_record[1] in ['M','=']:
        new_read_rec=read_rec-start_dis
        new_start_dis=0
        return [new_read_rec,new_start_dis]
    else:
        return [read_rec,start_dis]

def chop_pacbio_read_by_pos(bam_in_new,chrom,start,end,flank_length):
    fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,chrom,start,end))
    out=[]
    for line in fbam:
        pbam=line.strip().split()
        if not pbam[0]=='@': 
            if int(pbam[3])<start+1:
                align_info=cigar2alignstart_by_pos(pbam[5],int(pbam[3]),start,end)
                align_start=align_info[0]
                miss_bp=align_info[1]
                if not miss_bp>flank_length/2:
                    target_read=pbam[9][align_start:]
                    if len(target_read)>end-start-miss_bp:                    
                        out.append([target_read[:end-start-miss_bp],miss_bp,pbam[0]])
    fbam.close()
    return out

def chromos_readin(ref):
    fin=open(ref+'.fai')
    chromos=[]
    for line in fin:
            pin=line.strip().split()
            chromos.append(pin[0])
    fin.close()
    return chromos

def chr_start_end_extract(pin):
    out=[pin[0],int(pin[1])]
    for x in pin[7].split(';'):
        if x[:4]=='END=' and x.split('=')[0]=='END':
            out.append(int(x.split('=')[1]))
    return out

def cluster_range_decide(other_cluster):
    out=[]
    for x in other_cluster:
        out.append([])
        out[-1].append([min(x[0]),max(x[0])])
        out[-1].append([min(x[1]),max(x[1])])
    return out

def cluster_size_decide(range_cluster):
    out=[]
    for x in range_cluster:
        size=(x[0][1]-x[0][0])*(x[1][1]-x[1][0])
        out.append(np.sqrt(size))
    return out

def clean_dotdata_m1(ref_dotdata):
    #cluster dots based on their distance to diagnal
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]-x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot_cluster=dis_cluster([i[0] for i in temp],40)
        for y in kept_dot_cluster:
            kept_dot.append([temp[i] for i in y])
    kept_region=[[x[0],x[1]] for x in kept_dot]
    out=[]
    for x in kept_dot:    out+=x
    return [out,kept_region]

def clean_dotdata_diagnal_m1b(ref_dotdata):
    #cluster dots based on their distance to diagnal
    if ref_dotdata==[]:    return [[],[]]
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]-x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot+=temp
    kept_region=[[x[0],x[1]] for x in kept_dot]
    return [kept_dot,kept_region]

def clean_dotdata_anti_diagnal_m1b(ref_dotdata):
    #cluster dots based on their distance to diagnal
    if ref_dotdata==[]:    return [[],[]]
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]+x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot+=temp
    kept_region=[[x[0],x[1]] for x in kept_dot]
    return [kept_dot,kept_region]

def clean_dotdata_diagnal_and_anti_diagnal(ref_dotdata):
    if ref_dotdata==[]:    return [[],[]]
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]-x_axis[i] for i in range(len(x_axis))]
    dis_to_anti_diagnal=[y_axis[i]+x_axis[i] for i in range(len(x_axis))]
    #[kept_dot_x,removed_dot_x]=dis_cluster_2(x_axis,dis_cff=10)
    #[kept_dot_y,removed_dot_y]=dis_cluster_2(y_axis,dis_cff=10)
    [kept_dot_diag,removed_dot_diag]=dis_cluster_2(dis_to_diagnal,dis_cff=10)
    [kept_dot_anti_diag,removed_dot_anti_diag]=dis_cluster_2(dis_to_anti_diagnal,dis_cff=10)
    rec=-1
    kept_dot=[]
    for x in ref_dotdata:
        rec+=1
        if rec in removed_dot_diag and rec in removed_dot_anti_diag:   continue
        else:   kept_dot.append(x)
    return kept_dot

def clean_dotdata_x_and_y(list_dotdata):
    [kept_dot_x,removed_dot_x]=dis_cluster_2([i[0] for i in list_dotdata],dis_cff=10)
    kept_list=[]
    for i in kept_dot_x:
        kept_list+=[list_dotdata[j] for j in i]
    [kept_dot_y,removed_dot_y]=dis_cluster_2([i[1] for i in kept_list],dis_cff=10)
    out_list=[]
    for i in kept_dot_y:
        out_list+=[list_dotdata[j] for j in i]
    return out_list

def clean_dotdata_m2(ref_dotdata):
    #keep dots closest to diagnal
    data_hash={}
    for x in ref_dotdata:
        if not x[0] in data_hash:    data_hash[x[0]]=x[1]
        else:
            if abs(x[1]-x[0])<abs(data_hash[x[0]]-x[0]):    data_hash[x[0]]=x[1]
    out=[[[i,data_hash[i]] for i in sorted(data_hash.keys())],[]]
    return out

def complementary(seq):
    seq2=[]
    for i in seq:
            if i in 'ATGCN':
                    seq2.append('ATGCN'['TACGN'.index(i)])
            elif i in 'atgcn':
                    seq2.append('atgcn'['tacgn'.index(i)])
    return ''.join(seq2)

def compute_bic(kmeans,X):
    """
    Computes the BIC metric for a given clusters
    Parameters:
    -----------------------------------------
    kmeans:  List of clustering object from scikit learn
    X     :  multidimension np array of data points
    Returns:
    -----------------------------------------
    BIC value
    """
    # assign centers and labels
    centers = [kmeans.cluster_centers_]
    labels  = kmeans.labels_
    #number of clusters
    m = kmeans.n_clusters
    # size of the clusters
    n = np.bincount(labels)
    #size of data set
    N, d = X.shape
    #compute variance for all clusters beforehand
    cl_var=[]
    for i in range(m):
        if not n[i] - m==0:
            cl_var.append((1.0 / (n[i] - m)) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
        else:
            cl_var.append(float(10**20) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
    const_term = 0.5 * m * calcu_log10(N)
    removed_indices = find_removed_indices_with_negative(cl_var)
    #print(n, N, d, const_term, cl_var, removed_indices)
    n =  [arr for i, arr in enumerate(n) if i not in removed_indices]
    cl_var = [arr for i, arr in enumerate(cl_var) if i not in removed_indices]
    BIC = np.sum([n[i] * calcu_log10(n[i]) -
           n[i] * calcu_log10(N) -
         ((n[i] * d) / 2) * calcu_log10(2*np.pi) -
          (n[i] / 2) * calcu_log10(cl_var[i]) -
         ((n[i] - m) / 2) for i in range(len(n))]) - const_term
    return(BIC)

def find_removed_indices_with_negative(arrays):
    removed_indices = []
    for i, arr in enumerate(arrays):
        arrays[i] = [0.0 if x == -0.0 else x for x in arr]
        if any(x < 0 for x in arrays[i]):
            removed_indices.append(i)
    return removed_indices

def dup_inv_ref_alt_bps_produce(sv_info,flank_length,alt_structure):
    bp_info=sorted(sv_info[1:3]+[sv_info[4]])
    block_len=bp_to_block_len([sv_info[0]]+bp_info)
    ref_bps=[bp_info[0]-flank_length]+bp_info+[bp_info[-1]+flank_length]
    alt_bps=ref_bps[:2]
    for let_single in alt_structure:
        alt_bps.append(alt_bps[-1]+block_len[let_single[0]])
    alt_bps+=[alt_bps[-1]+flank_length]
    return [ref_bps,alt_bps]

def dup_inv_dup_bps_produce(sv_info,flank_length,alt_structure):
    [ref_bps,alt_bps]=dup_inv_ref_alt_bps_produce(sv_info,flank_length,alt_structure)
    alt_bps_new=[i-alt_bps[0] for i in alt_bps]
    if len(alt_structure)==2:
        return [alt_bps_new[1:3],alt_bps_new[2:4]]
    else:
        return [alt_bps_new[1:3],alt_bps_new[3:5]]

def dotdata(kmerlen,seq1, seq2):
    nth_base = 1
    inversions = True
    hits = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
    return hits

def dis_cluster(dis_to_diagnal,dis_cff=10):
    #eg of dis_list=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 0, 683, 0, 683, 0, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 1025, 0, 341, 512, 683, 1025, 0, 170, 341, 512, 683, 1025, 0, 170, 341, 683, 854, 1025, 0, 170, 683, 854, 1025, 0, 683, 1025, 0, 0, 0, 0, 0, 0, 0, 0, 512, 0, 512, 0, 512, 1196, 0, 170, 341, 512, 854, 1196, 1367, 0, 170, 341, 512, 854, 1196, 1367, 0]
    sorted_list=sorted(dis_to_diagnal)
    sub_group1=[[sorted_list[0]]]
    for i in sorted_list[1:]:
        if i-sub_group1[-1][-1]<dis_cff:
            sub_group1[-1].append(i)
        else:
            sub_group1.append([i])
    remove_noise_1=[i for i in sub_group1 if len(i)>50]
    if remove_noise_1==[]:
        lenght_list=[len(i) for i in sub_group1]
        remove_noise_1=[i for i in sub_group1 if len(i)==max(lenght_list)]
    return [[i for i in range(len(dis_to_diagnal)) if dis_to_diagnal[i] in j] for j in remove_noise_1]

def dis_cluster_2(dis_to_diagnal,dis_cff=10):
    #eg of dis_list=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 0, 683, 0, 683, 0, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 1025, 0, 341, 512, 683, 1025, 0, 170, 341, 512, 683, 1025, 0, 170, 341, 683, 854, 1025, 0, 170, 683, 854, 1025, 0, 683, 1025, 0, 0, 0, 0, 0, 0, 0, 0, 512, 0, 512, 0, 512, 1196, 0, 170, 341, 512, 854, 1196, 1367, 0, 170, 341, 512, 854, 1196, 1367, 0]
    sorted_list=sorted(dis_to_diagnal)
    sub_group1=[[sorted_list[0]]]
    for i in sorted_list[1:]:
        if i-sub_group1[-1][-1]<dis_cff:
            sub_group1[-1].append(i)
        else:
            sub_group1.append([i])
    remove_noise_1=[i for i in sub_group1 if len(i)>10]
    out1=[[i for i in range(len(dis_to_diagnal)) if dis_to_diagnal[i] in j] for j in remove_noise_1]
    out_total=[]
    for i in out1:  out_total+=i
    out2=[i for i in range(len(dis_to_diagnal)) if not i in out_total]
    return [out1,out2]

def dis_to_diagnal_most_abundant_defined(alt_clean_dotdata):
    dis_to_diagnal_list=[i[1]-i[0] for i in alt_clean_dotdata]
    dis_range=[min(dis_to_diagnal_list)+i*float(max(dis_to_diagnal_list)-min(dis_to_diagnal_list))/10.0 for i in range(11)]
    kept_list1=find_longest_list(number_cluster(dis_to_diagnal_list,dis_range))
    kept_list2=[]
    for km in kept_list1:
        j=[min(km)+i*float(max(km)-min(km))/10.0 for i in range(11)]
        kept_list2+=find_longest_list(number_cluster(km,j))
    if len(kept_list2)==1:  return np.median(kept_list2[0])
    else:                   return 0

def dot_to_line(list_dotdata,gap=50,length=10):
    ####This function try to recognize lines in the dot plot
    dis_to_diagnal=[i[1]-i[0] for i in list_dotdata]
    cluster1=one_dimention_cluster_by_gap(dis_to_diagnal,gap,length)
    dot_cluster1=[[list_dotdata[i] for i in j] for j in cluster1]
    dis_to_anti_diagnal=[[i[1]+i[0] for i in j] for j in dot_cluster1]
    cluster2=[one_dimention_cluster_by_gap(i,gap,length) for i in dis_to_anti_diagnal]
    dot_cluster2=[ [[dot_cluster1[k][i] for i in j] for j in cluster2[k]] for k in range(len(cluster2))]
    out=[]
    for i in dot_cluster2:
        for j in i: out.append(j)
    return [[i[0],i[-1]] for i in out]

def dup_block_combine(dup_block,k1_hap,k2_hap):
    #eg of dup_block=['a', 'b'];     k1_hap='abcd'   ; k2_hap='abab'
    all_combines=[]
    for x in range(len(dup_block)):
        all_combines+=[''.join(list(i)) for i in list(itertools.combinations(dup_block,x+1))]
    all_combines=dup_block_combined_qc(all_combines)
    kept_dup=[]
    for x in all_combines[::-1]:
        if k2_hap.count(x)>1:
            kept_dup.append(x)
    return dup_block_kept_qc(kept_dup)[::-1]

def dup_block_combined_qc(all_combines):
    #eg of all_combines=['a', 'b', 'c', 'd', 'ab', 'ac', 'ad', 'bc', 'bd', 'cd', 'abc', 'abd', 'acd', 'bcd', 'abcd']
    out=[]
    for x in all_combines:
        if len(x)==1:   out.append(x)
        else:
            temp=[ord(i) for i in x]
            if interval_dis_calcu_max(temp)>1: continue
            else:   out.append(x)
    return out

def dup_block_kept_qc(kept_dup):
    #eg of kept_dup:
    out=[]
    if len(kept_dup)>0:
        out.append(kept_dup[0])
        for y in kept_dup[1:]:
            flag_y=0 
            for z in out:
                if y in z:  flag_y+=1
            if flag_y==0:   out.append(y)
    return out

def dup_let_recombind(vec_in):
    if vec_in==[]:
        return []
    else:
        vec2=sorted(vec_in)
        vec=[[vec2[0]]]
        for ka in vec2[1:]:
            if ord(ka)-ord(vec[-1][-1])==1:
                vec[-1].append(ka)
            else:
                vec.append([ka])
        vec3=[]
        for ka in vec:
            if len(ka)==1:
                vec3.append(ka)
            else:
                for kb in range(2,len(ka)+1):
                    for kc in ka[:(1-kb)]:
                        vec3.append([])
                        for kd in range(kb):
                            vec3[-1].append(ka[ka.index(kc)+kd])                
        vec4=[''.join(i) for i in vec3]
        return vec4

def editDistance(seq1,seq2,cost_matrix,move_matrix,r,c):
    [iCost, dCost , mCost]=[1,1,1]
    if cost_matrix[r][c]==float('inf'):
        if r==0 and c==0:
            cost_matrix[r][c]=0
        elif r==0:
            move_matrix[r][c]=0
            cost_matrix[r][c]=editDistance(seq1,seq2,cost_matrix,move_matrix,r,c-1)+iCost
        elif c==0:
            move_matrix[r][c]=1
            cost_matrix[r][c] = editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c)+dCost
        else:
            iDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r,c-1) + iCost
            dDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c) + dCost
            if seq1[r-1] == seq2[c-1]:
                mDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c-1)
            else:
                mDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c-1)+mCost
            if iDist<dDist:
                if iDist < mDist:  # insertion is optimal
                    move_matrix[r][c] = 0
                    cost_matrix[r][c] = iDist
                else:    #match is optimal
                    move_matrix[r][c] = 2
                    cost_matrix[r][c] = mDist
            else:
                if dDist < mDist:    #deletion is optimal
                    move_matrix[r][c]=1
                    cost_matrix[r][c]=dDist
                else:    #match is optimal
                    move_matrix[r][c]=2
                    cost_matrix[r][c]=mDist
    return cost_matrix[r][c]

def edit_dis_setup(seq1,seq2):
    cost_matrix=np.full((len(seq1)+1,len(seq2)+1),np.inf)
    move_matrix=np.full((len(seq1)+1,len(seq2)+1),-1.0)
    opt_dist=editDistance(seq1, seq2, cost_matrix,move_matrix, len(seq1), len(seq2))
    return opt_dist

def eu_dis_abs_calcu(ref_clean_dotdata):
    #eg of ref_clean_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    eu_dis_abs=[abs(i[0]-i[1]) for i in ref_clean_dotdata]
    return np.mean(eu_dis_abs)

def eu_dis_single_dot(dot_pos):
    #############################################
    ####eg of dot_pos=(1212, 1217)
    ####this function returns %deviation of the dot from diagnal
    #############################################
    if dot_pos[0]==0:   return abs(float(dot_pos[0]-dot_pos[1])/float(dot_pos[0]+1))
    else:               return abs(float(dot_pos[0]-dot_pos[1])/float(dot_pos[0]))

def eu_dis_dir_calcu(list_dotdata):
    #eg of list_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    eu_dis_abs=[i[0]-i[1] for i in list_dotdata if eu_dis_single_dot(i)>0.1]
    if eu_dis_abs==[]:  return 0.0001
    else:               return np.mean(eu_dis_abs)

def eu_dis_reg_calcu(list_dotdata):
    ratio_new=eu_y_vs_x_ratio_calcu(list_dotdata)
    eu_dis_abs=[ratio_new*i[0]-i[1] for i in list_dotdata if eu_dis_single_dot([ratio_new*i[0],i[1]])>0.15]
    if eu_dis_abs==[]:  return 0.0001
    else:               return abs(np.mean(eu_dis_abs))

def eu_dis_dots_within_10perc(ref_clean_dotdata):
    #eg of ref_clean_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    dis_10perc=[abs(float(i[0]-i[1])/float(i[0])) for i in ref_clean_dotdata if i[0]>0]
    return len([i for i in dis_10perc if i <0.16])    

def eu_dis_region_calcu(list_dotdata,list_bps):
    #eg of list_bps=[127906399, 127906515, 127907283, 127907399, 127907515]
    ref_bps_new=[i-list_bps[0] for i in list_bps]
    ref_region=[[] for i in range(len(ref_bps_new)-1)]
    reca=0
    recb=0
    while True:
        if reca==len(list_dotdata) or recb==len(ref_region): break
        if list_dotdata[reca][0]<ref_bps_new[recb+1]:            
            ref_region[recb].append(list_dotdata[reca])
            reca+=1
        else:
            recb+=1
    if reca<len(list_dotdata):
        ref_region[-1]+=list_dotdata[reca:]
    out=[eu_dis_dir_calcu(i) for i in ref_region]
    print (out)
    out_new=[i for i in out if abs(i)>1]
    if out_new==[]:     return 0.0001
    else:               return np.mean(out_new)

def eu_dis_reg_dup_block_calcu(list_dotdata,dup_block_bps):
    ref_region=[[] for i in range(len(dup_block_bps)+1)]
    for x in list_dotdata:
        if not x[0]<dup_block_bps[0][0] and not x[0]>dup_block_bps[0][1]:   ref_region[0].append(x)
        elif not x[0]<dup_block_bps[1][0] and not x[0]>dup_block_bps[1][1]: ref_region[1].append(x)
        else:                                                               ref_region[2].append(x)
    out=[eu_dis_dir_calcu(i) for i in ref_region]
    out[-1]=abs(out[-1])
    out_new=[i for i in out if abs(i)>1]
    if out_new==[]:     return 0.0001
    else:               return np.mean(out_new)

def eu_y_vs_x_ratio_single_dot(dot_pos):
    #############################################
    ####eg of dot_pos=(1212, 1217)
    ####this function returns %deviation of the dot from diagnal
    #############################################
    if dot_pos[0]==0:   return 1
    else:               return abs(float(dot_pos[1])/float(dot_pos[0]))

def eu_y_vs_x_ratio_calcu(ref_clean_dotdata):
    y_x_ratio_list=[round(eu_y_vs_x_ratio_single_dot(i),2) for i in ref_clean_dotdata if eu_dis_single_dot(i)<0.15]
    if y_x_ratio_list==[]:  return 1
    else:
        if len(unify_list(y_x_ratio_list))>1:
            nparam_density = scipy.stats.gaussian_kde(y_x_ratio_list)
            test=scipy.optimize.fmin(lambda x:-nparam_density.pdf(x),1)
            if abs(test[0]-1) <0.15:        return test[0]
            else:   return 1
        else:
            return unify_list(y_x_ratio_list)[0]

def find_longest_list(list_set):
    length=[len(i) for i in list_set]
    out_pos=[i for i in range(len(list_set)) if len(list_set[i])==max(length)]
    out=[list_set[i] for i in out_pos]
    return unify_list(out)

def flank_length_calculate(bps):
    if int(bps[-1])-int(bps[1])<100:
        flank_length=(int(bps[-1])-int(bps[1]))
    else:
        if int(bps[-1])-int(bps[1])<500:
            flank_length=int(bps[-1])-int(bps[1])
        else:
            flank_length=500
    return flank_length

def genoCN_extract(pin):
    out=[0,0]
    rec_pos=-1
    if 'CN' in pin[8]:
        for x in pin[8].split(':'):
            rec_pos+=1
            if x=='CN': break
    geno=[i.split(':')[rec_pos] for i in pin[9:]]
    out=[0 if i=='2' else 1 for i in geno]
    return out

def genotype_extract(pin):
    out=[0,0]
    rec_pos=-1
    if 'GT' in pin[8]:
        for x in pin[8].split(':'):
            rec_pos+=1
            if x=='GT': break
    geno=[i.split(':')[rec_pos] for i in pin[9:]]
    for i in geno:
        if '/' in i:
            if i=='./.': out.append(1)
            else:   out.append(sum([int(j) for j in i.split('/')]))
        elif '|' in i:
            if i=='.|.': out.append(1)
            else:   out.append(sum([int(j) for j in i.split('|')]))
        elif i=='.':    out.append(1)
    return out

def INS_length_detect(pin):
    out=0
    for x in pin[7].split(';'):
        if 'SVLEN=' in x:
            out=int(x.split('=')[1])
    return out

def intersect(a, b):
    return ''.join(sorted(list(set(a) & set(b))))

def interval_dis_calcu_max(pos_check):
    #eg of pos_check=[97,98]
    if len(pos_check)>1:
        out=[pos_check[i+1]-pos_check[i] for i in range(len(pos_check)-1)]
        return max(out)
    else:
        return 'NA'

def kept_lines_size_filter(size_list,square_size=400):
    #eg of size_list=[(1582, 415), (1677, 506)]
    if abs((size_list[1][0]-size_list[0][0])*(size_list[1][1]-size_list[0][1]))>square_size:    return 'TRUE'
    else:                                                                                       return 'FALSE'

def k_means_cluster(data_list):
    if max(data_list[0])-min(data_list[0])>10 and max(data_list[1])-min(data_list[1])>10:
        array_diagnal=np.array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
        ks = list(range(1,min([5,len(data_list[0])+1])))
        KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
        KMeans_predict=[cluster.KMeans(n_clusters = i, init="k-means++").fit_predict(array_diagnal) for i in ks]
        BIC=[]
        BIC_rec=[]
        for x in ks:
            if KMeans_predict[x-1].max()<x-1: continue
            else:
                BIC_i=compute_bic(KMeans[x-1],array_diagnal)
                if abs(BIC_i)<10**8:
                    BIC.append(BIC_i)
                    BIC_rec.append(x)
        #BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
        #ks_picked=ks[BIC.index(max(BIC))]
        ks_picked=BIC_rec[BIC.index(max(BIC))]
        if ks_picked==1:
            return [data_list]
        else:
            out=[]
            std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
            whitened = whiten(array_diagnal)
            centroids, distortion=kmeans(whitened,ks_picked)
            idx,_= vq(whitened,centroids)
            for x in range(ks_picked):
                group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
                out.append(group1)
            return out
    else:
        return [data_list]

def k_means_cluster_Predict(data_list,info):
    array_diagnal=np.array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
    ks = list(range(1,len(info)))
    KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
    BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
    ks_picked=ks[BIC.index(max(BIC))]
    if ks_picked==1:
        return [data_list]
    else:
        out=[]
        std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
        whitened = whiten(array_diagnal)
        centroids, distortion=kmeans(whitened,ks_picked)
        idx,_= vq(whitened,centroids)
        for x in range(ks_picked):
            group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
            out.append(group1)
        return out

def key_modify(key):
    if 'R' in key:
            key=key.replace('R','N')
    if 'r' in key:
            key=key.replace('r','n')
    if 'Y' in key:
            key=key.replace('Y','N')
    if 'y' in key:
            key=key.replace('y','n')
    if 'S' in key:
            key=key.replace('S','N')
    if 's' in key:
            key=key.replace('s','n')
    if 'W' in key:
            key=key.replace('W','N')
    if 'w' in key:
            key=key.replace('w','n')
    if 'K' in key:
            key=key.replace('K','N')
    if 'k' in key:
            key=key.replace('k','n')
    if 'M' in key:
            key=key.replace('M','N')
    if 'm' in key:
            key=key.replace('m','n')
    if 'B' in key:
            key=key.replace('B','N')
    if 'b' in key:
            key=key.replace('b','n')
    if 'D' in key:
            key=key.replace('D','N')
    if 'd' in key:
            key=key.replace('d','n')
    if 'H' in key:
            key=key.replace('H','N')
    if 'h' in key:
            key=key.replace('h','n')
    if 'V' in key:
            key=key.replace('V','N')
    if 'v' in key:
            key=key.replace('v','n')
    return key

def kmerhits(seq1, seq2, kmerlen, nth_base=1, inversions=False):
    # hash table for finding hits
    lookup = {}
    # store sequence hashes in hash table
    #print ("hashing seq1...")
    seq1len = len(seq1)
    for i in range(seq1len - kmerlen + 1):
        key = seq1[i:i+kmerlen]
        for subkey in subkeys(key, nth_base, inversions):
            lookup.setdefault(subkey, []).append(i)
    # match every nth base by look up hashes in hash table
    #print ("hashing seq2...")
    hits = []
    for i in range(len(seq2) - kmerlen + 1):
        key = seq2[i:i+kmerlen]
        # only need to specify inversions for one seq
        for subkey in subkeys(key, nth_base, False):
            subhits = []
            if kmerlen>40:
                print(('Window size:' +str(kmerlen)))
                for k1 in list(lookup.keys()):
                    if edit_dis_setup(k1,subkey)<int(kmerlen/10)+1:
                        subhits+=lookup[k1]
            else:
                subhits=lookup.get(subkey,[])
            if subhits != []:
                # store hits to hits list
                for hit in subhits:
                    hits.append((i, hit))
                # break out of loop to avoid doubly counting
                # exact matches
                break
    return hits

def let_to_block_info(let,let_hash,chromos):
    #eg of let='ab'; eg of let_hash={'a': ['chrY', '10818935', '10819073'], 'b': ['chrY', '10819073', '10926507'], '+': ['chrY', '10926507', '10927007'], '-': ['chrY', '10818435', 10818935]}
    out=[]
    for i in let:
        if not i=='^':
            out+=let_hash[i]
    return(block_modify(out,chromos))

def letter_subgroup(k2_hap):
    #eg of k2_hap='ac^b^'
    inverted_sv=[]
    for x in k2_hap:
        if not x=='^':  inverted_sv.append(x)
        else:   inverted_sv[-1]+='^'
    inverted_sv_2=[]
    for x in inverted_sv:
        if inverted_sv_2==[]: inverted_sv_2.append(x)
        else:
            if not '^' in inverted_sv_2[-1] and not '^' in x and ord(x)-ord(inverted_sv_2[-1][-1])==1:  inverted_sv_2[-1]+=x
            elif '^' in inverted_sv_2[-1] and '^'  in x and ord(x[0])-ord(inverted_sv_2[-1][-2])==-1:   inverted_sv_2[-1]+=x
            else:   inverted_sv_2.append(x)
    inverted_sv_3=[]
    for i in inverted_sv_2:
        if not '^' in i:    inverted_sv_3.append(i)
        else:
            inverted_sv_3.append(i.replace('^','')[::-1]+'^')
    return inverted_sv_3

def letter_split(let):
    #eg of let='c^ba'
    out=[]
    for x in let:
        if not x=='^':  out.append(x)
        else:   out[-1]+=x
    return out

def list_unify(list):
    out=[]
    for i in list:
        if not i in out:    out.append(i)
    return out

def quality(hits):
    """determines the quality of a list of hits"""
    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000
    goodhits = []
    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2
        if lower < hit[1] < upper:
            goodhits.append(hit)
    return goodhits

def makeDotplot_subfigure(hits, title,figure_pos):
    #"""generate a dotplot from a list of hits filename may end in the following file extensions: *.ps, *.png, *.jpg"""
    if len(hits)==0: 
        #print hits
        return ' '
    x, y = zip(* hits)
    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000
    hits2 = quality(hits)
    xlib_range=int(float(max(x))/float(10**(len(str(max(x)))-1)))+1
    if xlib_range<3:
        xlib=[(i+1)*10**(len(str(max(x)))-1) for i in range(xlib_range)]
        xlib_new=[xlib[0]/2]
        for xi in range(len(xlib)-1):
            xlib_new.append(xlib_new[0]*(2*(xi+1)+1))
        xlib+=xlib_new
        xlib.sort()
    elif xlib_range<5:
        xlib=[(i+1)*10**(len(str(max(x)))-1) for i in range(xlib_range)]
    else:
        xlib=[(i+1)*2*10**(len(str(max(x)))-1) for i in range(int(xlib_range/2+1)+1)]
    plt.subplot(figure_pos)
    plt.plot(x, y,'+',color='r')
    plt.xticks(xlib, [str(i) for i in xlib])
    plt.title(title)
    plt.grid(False)
    #print "%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits)))
    # create plot

def make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name):
    nth_base = 1
    inversions = True
    if not best_read_rec=='':
        if not best_read_rec==[]:
            dotdata_record=[dotdata(window_size,ref_seq,ref_seq),    dotdata(window_size,alt_seq,alt_seq),    dotdata(window_size,best_read_rec[0],ref_seq[best_read_rec[1]:]),    dotdata(window_size,best_read_rec[0],alt_seq[best_read_rec[1]:])]
            [hits_ref_ref,hits_alt_alt,hits_ref,hits_alt]=dotdata_record
            if not [] in dotdata_record: 
                if  len(out_figure_name.split('/')[-1])>150:
                    out_figure_name='/'.join(out_figure_name.split('/')[:-1])+'/'+out_figure_name.split('/')[-1][:140]+'.'+out_figure_name.split('.')[-1]
                fig=plt.figure(plt_li)
                makeDotplot_subfigure(hits_ref_ref,'ref vs. ref',221)
                makeDotplot_subfigure(hits_alt_alt,'alt vs. alt',222)
                makeDotplot_subfigure(hits_ref,'read vs. ref',223)
                makeDotplot_subfigure(hits_alt,'read vs. alt',224)
                plt.savefig(out_figure_name)
                #plt.show()
                plt.close(fig)

def minimize_pacbio_read_list(x,ideal_list_length=20):
    #eg of x=chop_pacbio_read_by_pos(bam_in_new,chrom,start,end,flank_length)
    if len(x)>ideal_list_length:
        out=[]
        temp_hash={}
        for y in x:
            if not y[1] in list(temp_hash.keys()):    temp_hash[y[1]]=[y]
            else:    temp_hash[y[1]]+=[y]
        for y in sorted(temp_hash.keys()):
            if len(out)<ideal_list_length:    out+=temp_hash[y]
        return out[:ideal_list_length]
    else: return x

def number_cluster(dis_to_diagnal_list,dis_range):
    dis_clu=[[] for i in dis_range]
    reca=0
    recb=1
    dis_to_diagnal_list.sort()
    while True:
        if reca==len(dis_to_diagnal_list) or recb==len(dis_range): break
        if dis_to_diagnal_list[reca]<dis_range[recb]:
            dis_clu[recb-1].append(dis_to_diagnal_list[reca])
            reca+=1
        else:
            recb+=1
    if reca<len(dis_to_diagnal_list):
        dis_clu[-1]+=dis_to_diagnal_list[reca:]
    return dis_clu

def one_dimention_cluster_by_gap(dim1,gap,length):
    out_hash={}
    for k1 in range(len(dim1)):
        if not dim1[k1] in list(out_hash.keys()): out_hash[dim1[k1]]=[]
        out_hash[dim1[k1]].append(k1)
    sorted_keys=sorted(out_hash.keys())
    out=[[sorted_keys[0]]]
    for k1 in sorted_keys[1:]:
        if k1-out[-1][-1]>gap:  out.append([k1])
        else:   out[-1].append(k1)
    out2=[]
    for k1 in out:
        out2.append([])
        for k2 in k1:
            out2[-1]+=out_hash[k2]
    out2=[i for i in out2 if len(i)>length]
    return out2

def path_mkdir(path):
    if not os.path.isdir(path):
        os.system(r'''mkdir %s'''%(path))

def path_modify(path):
    if not path[-1]=='/':
        path+='/'
    return path

def polarity_detect(pin):
    out='+'
    for x in pin[7].split(';'):
        if 'MEIINFO=' in x:
            out=x.split(',')[-1]
    return out

def qual_check_repetitive_region(dotdata_qual_check):
    #qual_check_R_2(dotdata_qual_check,txt_file)  #not necessary for validator
    diagnal=0
    other=[[],[]]
    for x in dotdata_qual_check:
        if x[0]==x[1]:
            diagnal+=1 
        else:
            if x[0]>x[1]:
                other[0].append(x[0])
                other[1].append(x[1])
    if len(dotdata_qual_check)>0 and float(len(other[0]))/float(len(dotdata_qual_check))>0.1 and float(len(other[0]))/float(len(dotdata_qual_check))<0.5:
        other_cluster=X_means_cluster_reformat(other)
        range_cluster=cluster_range_decide(other_cluster)
        size_cluster=cluster_size_decide(range_cluster)
    else:
        size_cluster=[0]
    return [float(diagnal)/float(len(dotdata_qual_check)),size_cluster]

def reverse(seq):
    return seq[::-1]

def ref_ref_deviate_lines_calcu(list_dotdata):
    ref_kept_deviate=[i for i in list_dotdata if eu_dis_single_dot(i) > 0 and i[1]>i[0]]
    ref_left_wings=dot_to_line(ref_kept_deviate)
    ref_left_new=[]
    for x in ref_left_wings:
        ref_left_new.append(x)
        ref_left_new.append([[i[1],i[0]] for i in x])
    out=[]
    for x in ref_left_new:
        if x[0][0]<x[1][0]: out.append(x)
        else:               out.append([x[1],x[0]])
    return [i for i in out if kept_lines_size_filter(i)=='TRUE']

def ref_ref_deviate_lines_describe(list_dotdata):
    ref_ref_lines=ref_ref_deviate_lines_calcu(list_dotdata)
    out=[]
    for line in ref_ref_lines:
        ratio=round(float(line[1][1]-line[0][1])/float(line[1][0]-line[0][0]),0)
        intercept=round(np.mean([line[1][1]-line[1][0],line[0][1]-line[0][0]]),0)
        out.append([ratio,intercept,line[0][0],line[1][0]])
    return out

def ref_deviate_lines_calcu(list_dotdata):
    ref_kept_deviate=[i for i in list_dotdata if eu_dis_single_dot(i) > 0.15]
    ref_left_wings=dot_to_line(ref_kept_deviate)
    return    [i for i in ref_left_wings if kept_lines_size_filter(i)=='TRUE']

def ref_seq_readin(ref,chrom,start,end,reverse_flag='FALSE'):
    #reverse=='TRUE': return rev-comp-seq   ; if not specified, default as 'FALSE'
    #else: return original seq
    fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(start),int(end)))
    fref.readline().strip().split()
    seq=''
    while True:
            pref=fref.readline().strip().split()
            if not pref: break
            seq+=pref[0]
    fref.close()
    if reverse_flag=='FALSE':
        return seq
    else:
        return reverse(complementary(seq))

def result_organize_ins(info_list):
    #eg of info_list=[key_event,vapor_score_event]=['chr2_82961201', [-9.228366096827557, -106.46718851834126, -667.0858781654538, -38.56838396416415, -64.87185751169045, -147.77261544769615, -28.29536680099185, -25.378519434143666, -17.23542013374081, -113.00564782332029, -64.53043553409316]]
    if len(info_list[1])>0:
        pos_values=[i for i in info_list[1] if float(i)>0]
        neg_values=[i for i in info_list[1] if not float(i)>0]
        geno_value=float(len(pos_values))/float(len(pos_values)+len(neg_values))
        if not pos_values==[]:
            qual_value=np.mean(pos_values)
        else:
            qual_value=0
        return [info_list[0]]+[qual_value,geno_value,','.join([str(round(float(i),2)) for i in info_list[1]])]
    else:
        return [info_list[0]]+['NA' for i in range(3)]

def simple_del_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='a/a'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_del_haploid_decide(k1_hap,x))
    return out

def simple_del_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='ab'  ;    eg of k2_hap='b'
    if k1_hap==k2_hap: return 'FALSE'   #no alt
    if k2_hap=='': return [i for i in k1_hap]
    if '^' in k2_hap:   return 'FALSE'  #check if inv included
    dup_test=[k2_hap.count(x) for x in k2_hap]
    if max(dup_test)>1:     return 'FALSE'  #check if dup included
    if len(k2_hap)==1 and len(k1_hap)>1:    return letter_subgroup(''.join([i for i in k1_hap if not i in k2_hap]))   #del
    pos_compare=[ord(k2_hap[i+1])-ord(k2_hap[i]) for i in range(len(k2_hap)-1)]
    if min(pos_compare)<1: return 'FALSE'
    return letter_subgroup(''.join([i for i in k1_hap if not i in k2_hap]))

def simple_inv_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='ab^/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_inv_haploid_decide(k1_hap,x))
    return out

def simple_inv_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='ab'  ;    eg of k2_hap='b^a^'
    if not '^' in k2_hap:   return 'FALSE'      #if not block inverted
    if len(k2_hap.replace('^',''))==1 and len(k1_hap)==1:   return [i for i in k1_hap]
    dup_test=[k2_hap.count(i) for i in k2_hap if not i=='^']
    if max(dup_test)>1: return 'FALSE'
    inverted_sv_new=letter_subgroup(k2_hap)
    if ''.join([i.replace('^','') for i in inverted_sv_new])==k1_hap: return [i[:-1] for i in inverted_sv_new if '^' in i]
    else:   return 'FALSE'

def simple_tandup_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='abb/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_tandup_haploid_decide(k1_hap,x))
    return out

def simple_tandup_haploid_decide(k1_hap,k2_hap):
    if '^' in k2_hap:   return 'FALSE'
    dup_count=[k2_hap.count(i) for i in k1_hap]
    if min(dup_count)<1 or max(dup_count)<2:    return 'FALSE'  #deletion structure inside
    out=[]
    temp1=[]
    for x in k2_hap:
        if temp1==[]:   temp1.append(x)
        elif ord(x)-ord(temp1[-1][-1])==1:  temp1[-1]+=x
        else:   temp1.append(x)
    overlap_portion=[]
    overlap_count=[]
    for x in temp1:
        if out==[]:
            out.append(x)
        else:
            overlap=intersect(out[-1],x)
            if not len(overlap) >len(out[-1]) and not len(overlap)>len(x):
                if out[-1][-len(overlap):]==x[:len(overlap)]:
                    out[-1]+=x[len(overlap):]
                    if not overlap in overlap_portion:
                        overlap_portion.append(overlap)
                        overlap_count.append(2)
                    else:
                        overlap_count[overlap_portion.index(overlap)]+=1
                else:
                    out.append(x)
            else:
                out.append(x)
    if ''.join(out)==k1_hap:
        return [overlap_portion,overlap_count]
    return 'FALSE'

def simple_disdup_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='bab/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_disdup_haploid_decide(k1_hap,x))
    return out

def simple_disdup_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='abcd'    ;   eg of k2_hap='babdcd'
    if not '^' in k2_hap:
        if simple_tandup_haploid_decide(k1_hap,k2_hap)=='FALSE': 
            dup_dis=letter_subgroup(k2_hap)
            overlap=[intersect(dup_dis[i],dup_dis[i+1]) for i in range(len(dup_dis)-1)]
            if len(list_unify(overlap))==len(overlap):
                dup_count=[k2_hap.count(i) for i in k1_hap]
                if not min(dup_count)<1 and not max(dup_count)<2:    #deletion structure inside
                    dup_block=[k1_hap[i] for i in range(len(dup_count)) if dup_count[i]>1]
                    dup_block_combined=dup_block_combine(dup_block,k1_hap,k2_hap)
                    dis_dup_check=[]
                    no_dup_block=[]
                    for x in k2_hap: 
                        if not x in dup_block:
                            no_dup_block.append(k2_hap.index(x))
                    for x in dup_block_combined:
                        dis_dup_check.append([])
                        for y in range(len(k2_hap)-len(x)+1):
                            if k2_hap[y:(y+len(x))]==x:
                                dis_dup_check[-1].append(y)
                    original_pos=[]
                    for x in itertools.product(*dis_dup_check):
                        x_modify_new=x_to_x_modify_new(x,dup_block_combined)
                        temp_structure=[k2_hap[i] for i in sorted(x_modify_new+no_dup_block)]
                        if ''.join(temp_structure)==k1_hap:
                            original_pos+=list(x)
                    if len(original_pos)>0:
                        insert_pos=[]
                        for i in dis_dup_check:
                            for j in i:
                                if not j in original_pos:
                                    insert_pos.append(j)
                        k2_hap_new=['-']+[i for i in k2_hap]+['+']
                        insert_block=[]
                        pos_rec=-1
                        for i in insert_pos:
                            pos_rec+=1
                            if len(dup_block_combined[pos_rec])==1:
                                insert_block.append([k2_hap_new[i],k2_hap_new[i+1],k2_hap_new[i+2]])
                            else:
                                insert_block.append([k2_hap_new[i]]+k2_hap_new[(i+1):(i+len(dup_block_combined[pos_rec])+2)])
                        #insert_block=[[k2_hap_new[i],k2_hap_new[i+1],k2_hap_new[i+2]] for i in insert_pos]
                        return [dup_block_combined,insert_block]
    return 'FALSE'

def simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length):
    #eg of info = ['a/a','/','chr1', 101553562, 101553905]
    #pick read around left breakpoints
    #block_length={}
    #for x in range(len(info)-4):
    #    block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
    bam_in_new_list=bam_in_decide(bam_in,sv_info)
    if bam_in_new_list=='': return [[],[],[]]
    x=[]
    for bam_in_new in bam_in_new_list:
        x+=chop_pacbio_read_by_pos(bam_in_new, sv_info[0],int(sv_info[1])-flank_length,int(sv_info[1])+flank_length,flank_length)
    x=minimize_pacbio_read_list(x)
    return x

def simple_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length):
    #eg of info = ['a/a','/','chr1', 101553562, 101553905]
    #eg of sv_info=['chr1', 101553562, 101553905]
    bam_in_new_list=bam_in_decide(bam_in,sv_info)
    if bam_in_new_list=='': return [[],[],[]]
    x=[]
    for bam_in_new in bam_in_new_list:
        x+=chop_pacbio_read_by_pos(bam_in_new, sv_info[0],int(sv_info[1])-flank_length,int(sv_info[-1])+flank_length,flank_length)
    x=minimize_pacbio_read_list(x)
    return x

def subkeys(key, nth_base, inversions):
    subkeys_info = []
    key=key_modify(key)
    keylen = len(key)
    # speed tip from:
    # http://wiki.python.org/moin/PythonSpeed/PerformanceTips#String_Concatenation
    if nth_base == 1:
        subkeys_info = [key]
    elif nth_base != 0:
        for k in range(nth_base):
            substr_list = [key[j] for j in range(keylen) if (j % nth_base == k)]
            subkeys_info.append("".join(substr_list))
    else:
        # nth_base = 0 is a special case for third base mismatches
        # for every codon, only include the first 2 bases in the hash
        subkeys_info = ["".join([key[i] for i in range(len(key)) if i % 3 != 2])]
    if inversions:
        for i in range(len(subkeys_info)):
            subkeys_info.append("".join([invert_base[c] for c in reversed(subkeys_info[i])]))
    return subkeys_info

def svtype_extract(pin):
    svtype=''
    for x in pin[7].split(';'):
        if 'SVTYPE' in x:
            svtype=x.split('=')[1]
    if svtype=='':
        svtype=pin[4].replace('<','').replace('>','')
    return svtype

def sv_len_extract(pin):
    svtype=''
    for x in pin[7].split(';'):
        if 'SVLEN' in x:
            svtype=x.split('=')[1]
    if svtype=='':
        svtype=0
    return svtype

def sv_seq_extract(pin):
    seq=''
    for x in pin[7].split(';'):
        if x[:4]=='SEQ=':
            seq=x.split('=')[1]
    return seq

def sv_insert_point_define(pin):
    svtype=[0,0]
    for x in pin[7].split(';'):
        if 'insert_point=' in x:
            svtype=x.split('=')[1].split(':')
    if svtype==[0,0]:
        sv_info=pin[:2]
    return svtype

def take_off_symmetric_dots(list_dotdata):
    left_part=[list_dotdata[i] for i in range(int(len(list_dotdata)/2))]
    right_part=[list_dotdata[i][::-1] for i in [len(list_dotdata)-1-i for i in range(int(len(list_dotdata)/2))]]
    left_new=[i for i in left_part if eu_dis_single_dot(i)>0.15]
    right_new=[i for i in right_part if eu_dis_single_dot(i)>0.15]
    sym_dots=[]
    for i in left_new:
        for j in right_new:
            if abs(i[0]-j[0])<6 and abs(i[1]-j[1])<6:
                sym_dots.append(i)
                sym_dots.append(j[::-1])
    out_dots=[i for i in list_dotdata if not i in sym_dots]
    return out_dots

def two_dimention_cluster_by_gap(dim1,dim2,gap,length):
    ####dim1 and dim2 are both one-dimention data list,of the same length;
    ####dots with distance >gap on any dimention would be clustered into different group
    ####return position of each dots in each subcluster
    clu_1=one_dimention_cluster_by_gap(dim1,gap,length)
    clu_2_sub=[[dim2[i] for i in j] for j in clu_1]
    out=[]
    for i in clu_2_sub:
        out+=one_dimention_cluster_by_gap(i,gap,length)
    return out

def unify_list(list):
    out=[]
    for x in list:
        if not x in out:
            out.append(x)
    return out

def vapor_CANNOT_CLASSIFY_VapoR(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['ab_ab', 'b_b^', 'chr7', '70955990', '70961199', '70973901']
    ref_sv=sv_info[0].split('_')
    alt_sv=list_unify([i for i in sv_info[1].split('_') if not i in ref_sv])
    chromos=chromos_readin(ref)
    bp_info=block_subsplot(sv_info[2:],chromos)
    flank_length=max([flank_length_calculate(i) for i in bp_info])
    vapor_score_list=[]
    run_flag=0
    if len(bp_info)==1: #inter-chromosome event not considered here
        if bp_info[0][-1]-bp_info[0][1]<default_max_sv_test:
            ref_seq=ref_seq_readin(ref,bp_info[0][0],bp_info[0][1]-flank_length,bp_info[0][-1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,bp_info[0],flank_length)
                bp_let_hash=bp_to_chr_hash(bp_info[0],chromos,flank_length)
                if len(all_reads)>num_reads_cff:
                    run_flag+=1
                    bp_let_seq={}
                    for i in list(bp_let_hash.keys()):
                        bp_let_seq[i]=ref_seq_readin(ref,bp_let_hash[i][0],int(bp_let_hash[i][1]),int(bp_let_hash[i][-1]))
                    for alt_allele in alt_sv:
                        alt_seq=ref_seq[:flank_length]
                        for i in letter_split(alt_allele):
                            if not '^' in i:    alt_seq+=bp_let_seq[i]
                            else:               alt_seq+=reverse(complementary(bp_let_seq[i[0]]))
                        alt_seq+=ref_seq[-flank_length:]
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            best_read_rec=''
                            for x in all_reads:
                                if not max([alt_allele.count(i) for i in alt_allele]+[0])>1:    vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                                else:                                                           vapor_single_read_score=calcu_vapor_single_read_score_directed_dis_m1b_redefine_diagnal(ref_seq,alt_seq,x,window_size)
                                if not 0 in vapor_single_read_score:
                                    vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                    if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                            make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,'.'.join(out_figure_name.split('.')[:-1]+[ref_sv[0]+'.vs.'+alt_allele,out_figure_name.split('.')[-1]]))
        if run_flag==0:
            for alt_allele in alt_sv:
                alt_juncs=block_around_check(alt_allele,ref_sv[0])
                bp_let_hash=bp_to_chr_hash(bp_info[0],chromos,flank_length)
                for alt_jun in alt_juncs:
                    if not '^' in alt_jun[0]:                        ref_seq_a=ref_seq_readin(ref,bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]-flank_length,bp_let_hash[alt_jun[0][0]][2]+flank_length)
                    else:                                            ref_seq_a=reverse(complementary(ref_seq_readin(ref,bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]-flank_length,bp_let_hash[alt_jun[0][0]][1]+flank_length)))
                    if not '^' in alt_jun[1]:                        ref_seq_b=ref_seq_readin(ref,bp_let_hash[alt_jun[1][0]][0],bp_let_hash[alt_jun[1][0]][1]-flank_length,bp_let_hash[alt_jun[1][0]][1]+flank_length)
                    else:                                            ref_seq_b=reverse(complementary(ref_seq_readin(ref,bp_let_hash[alt_jun[1][0]][0],bp_let_hash[alt_jun[1][0]][2]-flank_length,bp_let_hash[alt_jun[1][0]][2]+flank_length)))
                    [window_size,window_size_qc]=window_size_refine(ref_seq_a+ref_seq_b)
                    if not window_size=='Error':
                        alt_seq=ref_seq_a[-flank_length:]+ref_seq_b[:flank_length]
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            if not '^' in alt_jun[0]:                   all_reads_a=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]],flank_length)
                            else:                                       all_reads_a=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]],flank_length)
                            #if not '^' in alt_jun[1]:                   all_reads_b=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]],flank_length)
                            #else:                                       all_reads_b=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]],flank_length)
                            if len(all_reads_a)>0:
                                for x in all_reads_a:
                                    vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq_a,alt_seq,x,window_size)
                                    if not 0 in vapor_single_read_score:
                                        vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            #if len(all_reads_b)>0:
                            #    for x in all_reads_a:
                            #        vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq_a,alt_seq,x,window_size)
                            #        if not vapor_single_read_score==[0,0] and not vapor_single_read_score[0]==0:
                            #            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
    return vapor_score_list

def vapor_del_inv_Vapor(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    sv_block=[sv_info[0][0],sv_info[0][1],sv_info[-1][2]]
    flank_length=flank_length_calculate(sv_block)
    vapor_score_list=[]
    best_read_rec=''
    if sv_info[1][1]-sv_info[0][2]<100:
        if sv_block[2]-sv_block[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
            ref_seq=ref_seq_readin(ref,sv_block[0],sv_block[1]-flank_length,sv_block[2]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq[:flank_length]
                for x in sv_info:
                    if x[-1]=='del':continue
                    elif x[-1]=='inv':  alt_seq+=reverse(complementary(ref_seq_readin(ref,x[0],x[1],x[2])))
                alt_seq+=ref_seq[-flank_length:]
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    all_reads=simple_chop_pacbio_read_simple_short(bam_in,sv_block[:2]+[sv_block[1]+len(alt_seq)-2*flank_length],flank_length)
                    if len(all_reads)>num_reads_cff:
                        best_read_rec=''
                        for x in all_reads:
                            vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not 0 in vapor_single_read_score:
                                vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                        make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
                    else:
                        if len(sv_info)==2 and [i[-1] for i in sv_info]==['del','inv']:
                            vapor_score_list=vapor_long_del_inv(bam_in,ref,sv_info,out_figure_name)
        else:
            if len(sv_info)==2 and [i[-1] for i in sv_info]==['del','inv']:
                vapor_score_list=vapor_long_del_inv(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name)
    else:
        for sub_sv_info in sv_info:
            if 'del' in sub_sv_info:    vapor_score_list+=vapor_simple_del_Vapor(bam_in,ref,sub_sv_info[:-1],'.'.join(out_figure_name.split('.')[:-1])+'_'.join([str(i) for i in sub_sv_info])+'.'+out_figure_name.split('.')[-1])
            elif 'inv' in sub_sv_info:    vapor_score_list+=vapor_simple_inv_Vapor(bam_in,ref,sub_sv_info[:-1],'.'.join(out_figure_name.split('.')[:-1])+'_'.join([str(i) for i in sub_sv_info])+'.'+out_figure_name.split('.')[-1])
    return vapor_score_list

def vapor_dup_inv_VapoR(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['chr1', 114103333, 114103408, 'chr1', 114111746]
    sv_info[1:3]=[int(i) for i in sv_info[1:3]]
    dup_block=sv_info[:3]
    ins_point=[sv_info[3],int(sv_info[4])]
    flank_length=flank_length_calculate(dup_block)
    vapor_score_list=[]
    if sv_info[0]==sv_info[3]:  #dup inv on the same chromosome
        bp_info=sorted(sv_info[1:3]+[sv_info[4]])
        run_flag=0
        if sv_info[0]==sv_info[3] and max(bp_info)-min(bp_info)<default_max_sv_test: #try to test on the whole region
            ref_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,max(bp_info)+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                run_flag+=1
                ref_structure=[i for i in 'ab']
                if sv_info[4]>sv_info[2]:           alt_structure=['a','b','a^']
                elif sv_info[4]<sv_info[1]:         alt_structure=['b^','a','b']
                else:                               alt_structure=['a','a^']
                #[ref_bps,alt_bps]=dup_inv_ref_alt_bps_produce(sv_info,flank_length,alt_structure)
                #dup_block_bps=dup_inv_dup_bps_produce(sv_info,flank_length,alt_structure)
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,[sv_info[0]]+bp_info+[bp_info[-1]+sv_info[2]-sv_info[1]],flank_length)
                if len(all_reads)>num_reads_cff:
                    alt_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,min(bp_info))
                    a_seq=ref_seq_readin(ref,sv_info[0],bp_info[0],bp_info[1])
                    b_seq=ref_seq_readin(ref,sv_info[0],bp_info[1],bp_info[2])
                    for x in alt_structure:
                        if x=='a':  alt_seq+=a_seq
                        elif x=='a^':   alt_seq+=reverse(complementary(a_seq))
                        elif x=='b':    alt_seq+=b_seq
                        elif x=='b^':   alt_seq+=reverse(complementary(b_seq))
                    alt_seq+=ref_seq_readin(ref,sv_info[0],max(bp_info),max(bp_info)+flank_length)
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        best_read_rec=''
                        for x in all_reads:
                            vapor_single_read_score=calcu_vapor_single_read_score_directed_dis_m1b_redefine_diagnal(ref_seq,alt_seq,x,window_size)
                            if not 0 in vapor_single_read_score and not math.isnan(vapor_single_read_score[0]) and not math.isnan(vapor_single_read_score[1]):
                                vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                        make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
        if run_flag==0:
            if max(bp_info)-min(bp_info)<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
                    if len(all_reads)>num_reads_cff:
                        alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[2])))+ref_seq[-flank_length:]
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            best_read_rec=''
                            for x in all_reads:
                                vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                                if not 0 in vapor_single_read_score and not math.isnan(vapor_single_read_score[0]) and not math.isnan(vapor_single_read_score[1]):
                                    vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                    if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                            make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
            else:
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
                    if len(all_reads)>num_reads_cff:
                        alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,dup_block[0],dup_block[2]-flank_length,dup_block[2])))
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            best_read_rec=''
                            for x in all_reads:
                                vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                                if not 0 in vapor_single_read_score and not math.isnan(vapor_single_read_score[0]) and not math.isnan(vapor_single_read_score[1]):
                                    vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                    if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                            make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vapor_long_del_inv(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=[['chr19', 46275941, 46314150, 'del'], ['chr19', 46314150, 46314312, 'inv']]
    vapor_score_list=[]
    best_read_rec=''
    flank_length=500
    ref_seq=ref_seq_readin(ref,sv_info[0][0],sv_info[0][1]-flank_length,sv_info[1][1]+flank_length) 
    [window_size,window_size_qc]=window_size_refine(ref_seq)
    if not window_size=='Error':
        alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,sv_info[1][0],sv_info[1][2]-flank_length,sv_info[1][2])))
        [window_size,window_size_qc]=window_size_refine(alt_seq)
        if not window_size=='Error':
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info[0],flank_length)
            if len(all_reads)>num_reads_cff:
                best_read_rec=''
                for x in all_reads:
                    vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    if not 0 in vapor_single_read_score:
                        vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def write_test_data(file_out,ref_dotdata,alt_dotdata):
    fo=open(file_out+'.ref','w')
    for i in ref_dotdata:    print('\t'.join([str(j) for j in i]), file=fo)
    fo.close()
    fo=open(file_out+'.alt','w')
    for i in alt_dotdata:    print('\t'.join([str(j) for j in i]), file=fo)
    fo.close()

def vapor_simple_del_Vapor(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['chr1', 101553562, 101553905]
    flank_length=flank_length_calculate(sv_info)
    vapor_score_list=[]
    best_read_rec=''
    if sv_info[2]-sv_info[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
        all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
        if len(all_reads)>num_reads_cff:
            ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[2]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq[:flank_length]+ref_seq[-flank_length:]
                best_read_rec=''
                for x in all_reads:
                    vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                    vapor_single_read_score_2=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    #for deletions, both calcu should be done to make sure, as I do have the over-pos issue
                    if not 0 in vapor_single_read_score and not 0 in vapor_single_read_score_2:
                        vapor_score_list.append(min([1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]) , 1-float(vapor_single_read_score_2[1])/float(vapor_single_read_score_2[0])]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    elif not 0 in vapor_single_read_score:
                        vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    elif not 0 in vapor_single_read_score_2:
                        vapor_score_list.append(1-float(vapor_single_read_score_2[1])/float(vapor_single_read_score_2[0]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x                       
                make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    else:
        #change the way of read in ref seq
        all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
        if len(all_reads)>num_reads_cff:
            ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1])+ref_seq_readin(ref,sv_info[0],sv_info[2],sv_info[2]+flank_length)
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    best_read_rec=''
                    for x in all_reads:
                        vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vapor_simple_tandup_Vapor(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
   #vapor_simple_del_Vapor(bam_in,ref,x[:-2],out_path+sample_name+'.TANDUP.'+key_event+'.png')
    flank_length=flank_length_calculate(sv_info)
    vapor_score_list=[]
    best_read_rec=''
    if sv_info[2]-sv_info[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
        ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[2]+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            alt_seq=ref_seq[:flank_length]+ref_seq[flank_length:(-flank_length)]+ref_seq[flank_length:(-flank_length)]+ref_seq[-flank_length:]
            [window_size,window_size_qc]=window_size_refine(alt_seq)
            if not window_size=='Error':
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,sv_info[:2]+[sv_info[1]+2*(sv_info[2]-sv_info[1])],flank_length)
                if len(all_reads)>num_reads_cff:
                    best_read_rec=''
                    for x in all_reads:
                        vapor_single_read_score=calcu_vapor_single_read_score_directed_dis_m1b_redefine_diagnal(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
                    return vapor_score_list
    ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[2]-flank_length,sv_info[2]+flank_length)
    [window_size,window_size_qc]=window_size_refine(ref_seq)
    if not window_size=='Error':
        alt_seq=ref_seq_readin(ref,sv_info[0],sv_info[2]-flank_length,sv_info[2])+ref_seq_readin(ref,sv_info[0],sv_info[1],sv_info[1]+flank_length)
        [window_size,window_size_qc]=window_size_refine(alt_seq)
        if not window_size=='Error':
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,[sv_info[0],sv_info[2]],flank_length)
            if len(all_reads)>num_reads_cff:
                best_read_rec=''
                for x in all_reads:
                    vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    if not 0 in vapor_single_read_score:
                        vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vapor_simple_disdup_Vapor(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    sv_info[1:3]=[int(i) for i in sv_info[1:3]]
    dup_block=sv_info[:3]
    ins_point=[sv_info[3],int(sv_info[4])]
    flank_length=flank_length_calculate(dup_block)
    vapor_score_list=[]
    best_read_rec=''
    bp_info=sorted([int(i) for i in sv_info[1:3]+[sv_info[4]]])
    run_flag=0
    if sv_info[0]==sv_info[3] and max(bp_info)-min(bp_info)<default_max_sv_test: #try to test on the whole region
        ref_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,max(bp_info)+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            all_reads=simple_chop_pacbio_read_simple_short(bam_in,[sv_info[0]]+bp_info+[int(bp_info[-1])+sv_info[2]-sv_info[1]],flank_length)
            if len(all_reads)>num_reads_cff:
                run_flag+=1
                ref_structure=[i for i in 'ab']
                if sv_info[4]>sv_info[2]:   alt_structure=['a','b','a']
                elif sv_info[4]<sv_info[1]: alt_structure=['b','a','b']
                alt_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,min(bp_info))
                a_seq=ref_seq_readin(ref,sv_info[0],bp_info[0],bp_info[1])
                b_seq=ref_seq_readin(ref,sv_info[0],bp_info[1],bp_info[2])
                for x in alt_structure:
                    if x=='a':  alt_seq+=a_seq
                    elif x=='b':    alt_seq+=b_seq
                alt_seq+=ref_seq_readin(ref,sv_info[0],max(bp_info),max(bp_info)+flank_length)
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    best_read_rec=''
                    for x in all_reads:
                        vapor_single_read_score=calcu_vapor_single_read_score_directed_dis_m1b_redefine_diagnal(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    if run_flag==0:
        if max(bp_info)-min(bp_info)<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
            if len(all_reads)>num_reads_cff:
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[2])+ref_seq[-flank_length:]
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        best_read_rec=''
                        for x in all_reads:
                            vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not 0 in vapor_single_read_score:
                                vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                        make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
        else:
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
            if len(all_reads)>num_reads_cff:
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[1]+flank_length)
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        best_read_rec=''
                        for x in all_reads:
                            vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                            if not 0 in vapor_single_read_score:
                                vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                                if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                        make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vapor_simple_ins_Vapor(num_reads_cff,plt_li,bam_in,ref,ins_pos,ins_seq,out_figure_name,POLARITY):
    #eg of ins_pos='chr1_83144055'
    #eg of bam_in='/nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'
    #eg of ref='/scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    if POLARITY=='+':    ins_seq_2=ins_seq
    elif POLARITY=='-':   ins_seq_2=reverse(complementary(ins_seq))
    flank_length=default_flank_length if len(ins_seq)>default_flank_length else len(ins_seq)
    vapor_score_list=[]
    best_read_rec=''
    all_reads=simple_chop_pacbio_read_simple_short(bam_in,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]]+[int(ins_pos.split('_')[-1])+len(ins_seq)],flank_length)
    if len(all_reads)>num_reads_cff:
        if len(ins_seq)<5000:
            ref_seq=ref_seq_readin(ref,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]][0],int(ins_pos.split('_')[-1])-flank_length,int(ins_pos.split('_')[-1])+flank_length+len(ins_seq),reverse_flag='FALSE')
            [window_size,window_size_qc]=window_size_refine(ref_seq+ins_seq)
        else:
            ref_seq=ref_seq_readin(ref,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]][0],int(ins_pos.split('_')[-1])-flank_length,int(ins_pos.split('_')[-1])+flank_length,reverse_flag='FALSE')
            [window_size,window_size_qc]=window_size_refine(ref_seq)        
        if not window_size=='Error':
            alt_seq=ref_seq_readin(ref,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]][0],int(ins_pos.split('_')[-1])-flank_length,int(ins_pos.split('_')[-1]),reverse_flag='FALSE')+ins_seq_2+ref_seq_readin(ref,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]][0],int(ins_pos.split('_')[-1]),int(ins_pos.split('_')[-1])+flank_length,reverse_flag='FALSE')
            best_read_rec=''
            if len(all_reads)>num_reads_cff:
                for x in all_reads:
                    if float(x[0].count('N')+x[0].count('n'))/float(len(x[0]))<0.1:
                        vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
            else:
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]]+[int(ins_pos.split('_')[-1])],flank_length)
                for x in all_reads:
                    if float(x[0].count('N')+x[0].count('n'))/float(len(x[0]))<0.1:
                        vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
            if ins_seq_2.count('X')==len(ins_seq_2):        make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,ref_seq[2:flank_length],out_figure_name)
            else:            make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vapor_simple_inv_Vapor(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name):
    #vapor_simple_inv_Vapor(bam_in,ref,y,out_path+sample_name+'.INV.'+key_event+'.png')
    #eg of sv_info=['chr1', 101553562, 101553905]
    flank_length=flank_length_calculate(sv_info)
    vapor_score_list=[]
    best_read_rec=''
    if sv_info[2]-sv_info[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
        ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[2]+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq[flank_length:(-flank_length)]))+ref_seq[-flank_length:]
            [window_size,window_size_qc]=window_size_refine(alt_seq)
            if not window_size=='Error':
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
                if len(all_reads)>num_reads_cff:
                    best_read_rec=''
                    for x in all_reads:
                        vapor_single_read_score=calcu_vapor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                        if not 0 in vapor_single_read_score:
                            vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                            if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                    make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
                    return vapor_score_list
    ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1]+flank_length)
    [window_size,window_size_qc]=window_size_refine(ref_seq)
    if not window_size=='Error':
        alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,sv_info[0],sv_info[2]-flank_length,sv_info[2],'TRUE')
        [window_size,window_size_qc]=window_size_refine(alt_seq)
        if not window_size=='Error':
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
            if len(all_reads)>num_reads_cff:
                best_read_rec=''
                for x in all_reads:
                    vapor_single_read_score=calcu_vapor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    if not 0 in vapor_single_read_score:
                        vapor_score_list.append(1-float(vapor_single_read_score[1])/float(vapor_single_read_score[0]))
                        if vapor_score_list[-1]==max(vapor_score_list):best_read_rec=x
                make_event_figure_1(plt_li,vapor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return vapor_score_list

def vcf_rec_hash_modify(vcf_rec_hash):
    out={}
    for k1 in list(vcf_rec_hash.keys()):
        if not vcf_rec_hash[k1] in list(out.keys()):  out[vcf_rec_hash[k1]]=[]
        out[vcf_rec_hash[k1]].append(k1)
    return out

def vcf_vapor_modify(vcf_input,vcf_rec_hash_new):
    vapor_input=vcf_input+'.vapor'
    vapor_rec={}
    fin=open(vcf_input)
    vcf_info_hash={}
    rec=-1
    for line in fin:
        rec+=1
        pin=line.strip().split()
        if not pin[0][0]=='#':
            #if pin[6]=='PASS':
                vcf_info_hash[rec]=pin
    fin.close()
    fin=open(vapor_input)
    keep_rec=[]
    for line in fin:
        pin=line.strip().split()
        if pin[0] in list(vcf_rec_hash_new.keys()):   
            sv_rec_label=vcf_rec_hash_new[pin[0]]
            for y in sv_rec_label:
                vcf_info_hash[y]+=[round(float(i),2) if not i=='NA' else i for i in pin[1:3]]+[pin[3]]+[round(float(pin[4]),2) if not pin[4]=='NA' else pin[4]]+[pin[5]]
                keep_rec.append(y)
    fin.close()
    fo=open(vapor_input,'w')
    print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE','VaPoR_gs','VaPoR_GT','VaPoR_GQ','VaPoR_Rec']), file=fo)
    for k1 in sorted(vcf_info_hash.keys()):
        if k1 in keep_rec:
            print('\t'.join([str(i) for i in vcf_info_hash[k1]]), file=fo)
    fo.close()

def vcf_vapor_modify(vcf_input,vcf_rec_hash_new):
    vapor_input=vcf_input+'.vapor'
    vapor_rec={}
    fin=open(vcf_input)
    vcf_info_hash={}
    ## Store meta-information and header for VCF file
    meta_info=[]
    header=[]
    rec=-1
    for line in fin:
        pin=line.strip().split()
        if not pin[0][0]=='#':
            #if pin[6]=='PASS':
                rec+=1
                vcf_info_hash[rec]=pin
        ## Store meta-information and header VCF file
        elif not pin[0]=='#CHROM':
            meta_info.append(pin)
        else:
            header=pin
    fin.close()
    fin=open(vapor_input)
    keep_rec=[]
    for line in fin:
        pin=line.strip().split()
        if pin[0] in list(vcf_rec_hash_new.keys()):   
            sv_rec_label=vcf_rec_hash_new[pin[0]]
            for y in sv_rec_label:
                ##                GS=[round(float(i),2) if not i=='NA' else i for i in [pin[2]]]
                VaPoR_GS=round(float(pin[2]),2) if not pin[2]=='NA' else pin[2]
                VaPoR_GT=pin[3]
                VaPoR_GQ=round(float(pin[4]),2) if not pin[4]=='NA' else pin[4]
                VaPoR_REC=pin[5]
                ##                vcf_info_hash[y]+=[round(float(i),2) if not i=='NA' else i for i in pin[1:3]]+[pin[3]]+[round(float(pin[4]),2) if not pin[4]=='NA' else pin[4]]+[pin[5]]
                vcf_info_hash[y][7]+=';VaPor_GS='+str(VaPoR_GS)+';VaPor_GT='+str(VaPoR_GT)+';VaPor_GQ='+str(VaPoR_GQ)+';VaPor_REC='+str(VaPoR_REC)
                keep_rec.append(y)
    fin.close()
    fo=open(vapor_input,'w')
    prev_meta_data=''
    current_meta_data=''
    ## Add vcf meta-information to the file
    for line in meta_info:
        joined_line=' '.join(line)
        current_meta_data=joined_line.split('=')[0]
        if prev_meta_data=='##INFO' and not current_meta_data=='##INFO':
            ## Add additional meta-info to file
            print('##INFO=<ID=VaPoR_GS,Number=1,Type=Float,Description="VaPoR Score, representing the percentage of transverse long reads that support the prediction">', file=fo)
            print('##INFO=<ID=VaPoR_GT,Number=1,Type=String,Description="Genotype with the highest likelihood as estimated by VaPoR">', file=fo)
            print('##INFO=<ID=VaPoR_GQ,Number=1,Type=Float,Description="Genotype quality score - likelihood of the second most likely genotype on a -log10 normalized scale"', file=fo)
            print('##INFO=<ID=VaPoR_REC,Number=.,Type=Float,Description="Similarity scores assigned to each of the reads traversings the predicted SV">', file=fo)
        print(joined_line, file=fo)
        prev_meta_data=current_meta_data 
    print('\t'.join(header), file=fo)
    for k1 in sorted(vcf_info_hash.keys()):
        if k1 in keep_rec:
            print('\t'.join([str(i) for i in vcf_info_hash[k1]]), file=fo)
    fo.close()

def window_size_refine(seq2,region_QC_Cff=0.4):
    window_size=10
    seq2=''.join([i for i in seq2 if not i=='X'])
    if not seq2.count('N')+seq2.count('n')>100:
        dotdata_qual_check=dotdata(window_size,seq2,seq2)
        if len(dotdata_qual_check)>0:
            region_QC=qual_check_repetitive_region(dotdata_qual_check)
            while True:
                if window_size>30: break
                if region_QC[0]>region_QC_Cff or sum(region_QC[1])/float(len(seq2))<0.3: break
                else:
                    window_size+=10
                    dotdata_qual_check=dotdata(window_size,seq2,seq2)
                    region_QC=qual_check_repetitive_region(dotdata_qual_check)
            return [window_size,region_QC]
        else:   return ['Error','Error']
    else:       return ['Error','Error']

def write_dotdata_to_file(file_out,dotdata_list):
    fo=open(file_out,'w')
    for k1 in dotdata_list:
        print(' '.join([str(i) for i in k1]), file=fo)
    fo.close()

def gt_estimate_log_likelihood(vapor_result):
    read_score_list=[float(i) for i in vapor_result[-1].split(',')]
    k=len(read_score_list)
    l=len([i for i in read_score_list if not i>0])
    m=2
    likelihood_0_0=log_likelihood_calcu(k,l,m,2)    #g is the number of ref alleles 
    likelihood_0_1=log_likelihood_calcu(k,l,m,1)
    likelihood_1_1=log_likelihood_calcu(k,l,m,0)
    gt_list=['0/0','0/1','1/1']
    gt_score=[likelihood_0_0,likelihood_0_1,likelihood_1_1]
    gt_ori_scale=[np.exp(i-max(gt_score)) for i in gt_score]
    gt_norm=[i/sum(gt_ori_scale) for i in gt_ori_scale]
    gt_qual=-np.log(np.median(gt_norm))/np.log(10)
    gt_out=gt_list[gt_score.index(max(gt_score))]
    if gt_out=='0/0' and vapor_result[-2]>.15 : gt_out='0/1'
    return [gt_out,gt_qual]

def log_likelihood_calcu(k,l,m,g,err=0.05):
    out=-k*np.log(m)
    for j in range(l):
        out+=np.log((m-g)*err + g*(1-err))
    for j in range(k-l):
        out+=np.log( (m-g)*(1-err) + g*err )
    return out 

def write_output_initiate(out_name):
    fo=open(out_name,'w')
    print('\t'.join(['#CHR','POS','END','SVTYPE','SVID','VaPoR_QS','VaPoR_GS','VaPoR_GT','VaPoR_GQ','VaPoR_Rec']), file=fo)
    fo.close()

def write_output_main(out_name,out_list):
    fo=open(out_name,'a')
    if not 'NA' in out_list:    print('\t'.join([str(i) for i in out_list[:-1]+gt_estimate_log_likelihood(out_list)+[out_list[-1]]]), file=fo)
    else:                       print('\t'.join([str(i) for i in out_list[:-1]+['NA','NA','NA']]), file=fo)
    fo.close()

def x_to_x_modify_new(x,dup_block_combined):
    x_modify=[[i] for i in list(x)]
    block_rec=-1
    for y in dup_block_combined:
        block_rec+=1
        if len(y)>1:
            x_modify[block_rec]+=[x_modify[block_rec][0]+1+i for i in range(len(y)-1)]
    x_modify_new=[]
    for y in x_modify:  x_modify_new+=y
    return x_modify_new

def X_means_cluster(data_list):
    temp_result=[i for i in k_means_cluster(data_list) if not i==[[],[]]]
    if temp_result==[data_list]:
        return temp_result[0]
    else:
        out=[]
        for i in temp_result:
            out+=X_means_cluster(i)
        return out

def X_means_cluster_reformat(data_list):
    out=X_means_cluster(data_list)
    out2=[]
    for y in range(int(len(out)/2)):
        out2.append([out[2*y],out[2*y+1]])
    return out2
